(function () {
  "use strict";

  const MAX_SAMPLE_READS = 10000;
  const PREFIX_BYTES = 64 * 1024;
  const MAX_READ_WINDOW = 1024 * 1024;
  const SAMPLE_BATCH_SIZE = 64;

  function setStatus(input, message, state) {
    input.dataset.nxStatus = state;
    const run = input.id.replace("sequencing_summary_file_", "");
    const status = document.getElementById(`sequencing_summary_status_${run}`);
    if (status) status.textContent = message;
  }

  async function calculateNx(file, input) {
    const decoder = new TextDecoder("utf-8");
    const lengthCounts = new Map();
    let totalBases = 0;
    let readCount = 0;
    const prefixEnd = Math.min(file.size, PREFIX_BYTES);
    const prefix = new Uint8Array(await file.slice(0, prefixEnd).arrayBuffer());
    const headerEnd = prefix.indexOf(10);
    if (headerEnd < 0) {
      throw new Error("Could not find a complete sequencing-summary header.");
    }
    const header = decoder.decode(prefix.subarray(0, headerEnd))
      .replace(/^\uFEFF/, "").replace(/\r$/, "").split("\t");
    const lengthColumn = header.indexOf("sequence_length_template");
    if (lengthColumn < 0) {
      throw new Error("Missing sequence_length_template column.");
    }

    const bodyStart = headerEnd + 1;
    const bodySize = file.size - bodyStart;
    if (bodySize <= 0) {
      throw new Error("The selected summary has no read rows.");
    }

    let prefixLineStart = bodyStart;
    let prefixLineCount = 0;
    let prefixLineBytes = 0;
    for (let position = bodyStart; position < prefix.length && prefixLineCount < 100; position++) {
      if (prefix[position] === 10) {
        prefixLineBytes += position + 1 - prefixLineStart;
        prefixLineStart = position + 1;
        prefixLineCount += 1;
      }
    }
    if (prefixLineCount === 0) {
      throw new Error("Could not estimate sequencing-summary row size.");
    }

    const averageLineBytes = prefixLineBytes / prefixLineCount;
    const estimatedReadCount = Math.max(1, Math.ceil(bodySize / averageLineBytes));
    const sampleCount = Math.min(MAX_SAMPLE_READS, estimatedReadCount);
    const initialReadWindow = Math.max(
      32,
      Math.min(64 * 1024, Math.ceil(averageLineBytes * 2))
    );
    const positions = Array.from({ length: sampleCount }, (_, index) =>
      bodyStart + Math.floor((index + 0.5) * bodySize / sampleCount)
    );
    const sampledLineStarts = new Set();

    async function readSampledLength(position, sampleIndex) {
      let windowSize = initialReadWindow;
      while (windowSize <= MAX_READ_WINDOW) {
        let start = Math.max(bodyStart, position - 1);
        const end = Math.min(file.size, position + windowSize);
        let bytes = new Uint8Array(await file.slice(start, end).arrayBuffer());
        let lineStart;
        if (position === bodyStart || bytes[position - start - 1] === 10) {
          lineStart = position;
        } else {
          const nextNewline = bytes.indexOf(10, position - start);
          if (nextNewline < 0) {
            if (end === file.size) {
              let backwardWindow = windowSize;
              while (true) {
                start = Math.max(bodyStart, position - backwardWindow);
                bytes = new Uint8Array(await file.slice(start, end).arrayBuffer());
                const previousNewline = bytes.lastIndexOf(10, position - start - 1);
                if (previousNewline >= 0) {
                  lineStart = start + previousNewline + 1;
                  break;
                }
                if (start === bodyStart) {
                  lineStart = bodyStart;
                  break;
                }
                backwardWindow *= 2;
                if (backwardWindow > MAX_READ_WINDOW) {
                  throw new Error("A sampled summary row exceeded the 1 MiB row-size limit.");
                }
              }
            } else {
              windowSize *= 2;
              continue;
            }
          } else {
            lineStart = start + nextNewline + 1;
          }
        }

        let lineEnd = bytes.indexOf(10, lineStart - start);
        if (lineEnd < 0 && end < file.size) {
          windowSize *= 2;
          continue;
        }
        if (sampledLineStarts.has(lineStart)) return null;
        sampledLineStarts.add(lineStart);

        const line = decoder.decode(bytes.subarray(
          lineStart - start,
          lineEnd < 0 ? bytes.length : lineEnd
        )).replace(/\r$/, "");
        const rawLength = (line.split("\t")[lengthColumn] || "").trim();
        if (!/^\d+$/.test(rawLength)) {
          throw new Error(`Invalid read length near sampled row ${sampleIndex + 1}.`);
        }
        const length = Number(rawLength);
        if (!Number.isSafeInteger(length)) {
          throw new Error(`Read length near sampled row ${sampleIndex + 1} is out of range.`);
        }
        return length;
      }
      throw new Error("A sampled summary row exceeded the 1 MiB row-size limit.");
    }

    for (let batchStart = 0; batchStart < positions.length; batchStart += SAMPLE_BATCH_SIZE) {
      const batch = positions.slice(batchStart, batchStart + SAMPLE_BATCH_SIZE);
      const sampledLengths = await Promise.all(
        batch.map((position, offset) => readSampledLength(position, batchStart + offset))
      );
      for (const length of sampledLengths) {
        if (!length) continue;
        lengthCounts.set(length, (lengthCounts.get(length) || 0) + 1);
        totalBases += length;
        if (!Number.isSafeInteger(totalBases)) {
          throw new Error("Sampled read bases exceed the browser's exact integer range.");
        }
        readCount += 1;
      }
      const percent = Math.floor(
        Math.min(batchStart + batch.length, positions.length) / positions.length * 100
      );
      setStatus(input, `Sampling locally: ${percent}%`, "processing");
      await new Promise(resolve => setTimeout(resolve, 0));
    }

    if (!readCount || !totalBases) {
      throw new Error("No positive read lengths were found in the sampled reads.");
    }

    const orderedLengths = Array.from(lengthCounts.keys()).sort((a, b) => b - a);
    const values = [];
    let cumulativeBases = 0;
    let percentile = 5;
    for (const length of orderedLengths) {
      cumulativeBases += length * lengthCounts.get(length);
      while (percentile <= 95 && cumulativeBases * 100 >= totalBases * percentile) {
        values.push({ x: percentile, y: Number((length / 1000).toFixed(3)) });
        percentile += 5;
      }
      if (percentile > 95) break;
    }
    if (values.length !== 19) {
      throw new Error("Could not calculate all Nx values.");
    }
    return { values, readCount, estimatedReadCount };
  }

  document.addEventListener("change", async function (event) {
    const input = event.target.closest("input.nx-summary-file");
    if (!input) return;

    // Prevent Shiny's built-in file binding from uploading the selected file.
    event.stopImmediatePropagation();
    event.preventDefault();

    const file = input.files && input.files[0];
    if (!file) return;
    const token = Number(input.dataset.nxToken || 0) + 1;
    input.dataset.nxToken = String(token);
    setStatus(input, "Reading locally…", "processing");

    try {
      const result = await calculateNx(file, input);
      if (Number(input.dataset.nxToken) !== token) return;
      const run = input.id.replace("sequencing_summary_file_", "");
      setStatus(
        input,
        `Ready: estimated Nx from ${result.readCount.toLocaleString()} sampled reads.`,
        "ready"
      );
      Shiny.setInputValue(
        `sequencing_summary_nx_${run}`,
        JSON.stringify(result.values),
        { priority: "event" }
      );
    } catch (error) {
      if (Number(input.dataset.nxToken) !== token) return;
      setStatus(input, error.message, "error");
      Shiny.setInputValue(
        input.id.replace("sequencing_summary_file_", "sequencing_summary_nx_"),
        null,
        { priority: "event" }
      );
    }
  }, true);

  ["dragenter", "dragover", "drop"].forEach(type => {
    document.addEventListener(type, function (event) {
      if (!event.target.closest(".nx-summary-control")) return;
      event.preventDefault();
      event.stopImmediatePropagation();
    }, true);
  });

  document.addEventListener("click", function (event) {
    if (event.target.closest("#reset")) {
      for (let run = 1; run <= 8; run += 1) {
        const input = document.getElementById(`sequencing_summary_file_${run}`);
        if (input) {
          input.dataset.nxToken = String(Number(input.dataset.nxToken || 0) + 1);
          input.value = "";
          setStatus(input, "Sampled locally; only estimated Nx values are sent to the server.", "");
        }
        Shiny.setInputValue(`sequencing_summary_nx_${run}`, null, { priority: "event" });
      }
      return;
    }

    if (!event.target.closest("#start")) return;
    const inputs = Array.from(document.querySelectorAll("input.nx-summary-file"));
    const states = inputs.map(input => input.dataset.nxStatus || "");
    const selected = states.filter(Boolean);
    if (selected.some(state => state !== "ready") ||
        (selected.length > 0 && selected.length !== inputs.length)) {
      event.preventDefault();
      event.stopImmediatePropagation();
      window.alert("Wait for each selected summary to finish, and process a summary for every run or none.");
    }
  }, true);

  Shiny.addCustomMessageHandler("resetNxSummary", function (run) {
    const input = document.getElementById(`sequencing_summary_file_${run}`);
    if (input) {
      input.dataset.nxToken = String(Number(input.dataset.nxToken || 0) + 1);
      input.value = "";
      setStatus(input, "Sampled locally; only estimated Nx values are sent to the server.", "");
    }
    Shiny.setInputValue(`sequencing_summary_nx_${run}`, null, { priority: "event" });
  });
})();

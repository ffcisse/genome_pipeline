(function () {
  "use strict";

  var DATA = DASHBOARD_DATA;
  var PRIMARY = DATA.meta.primary_grouping;
  var SUBGROUP = DATA.meta.subgroup_column;
  var GENOMES = DATA.genomes;
  var GENOME_NAMES = GENOMES.map(function (g) { return g.genome; });
  var PRIMARY_VALUES = DATA.grouping_values[PRIMARY];
  var SUBGROUP_VALUES = DATA.grouping_values[SUBGROUP];
  var PROPERTIES = DATA.meta.properties;
  var CDS_PROPERTIES = DATA.meta.cds_properties || [];

  var PALETTE_PRIMARY = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];
  var PALETTE_SUBGROUP = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"];

  var primaryColor = {};
  PRIMARY_VALUES.forEach(function (v, i) { primaryColor[v] = PALETTE_PRIMARY[i % PALETTE_PRIMARY.length]; });
  var subgroupColor = {};
  SUBGROUP_VALUES.forEach(function (v, i) { subgroupColor[v] = PALETTE_SUBGROUP[i % PALETTE_SUBGROUP.length]; });
  var genomeColor = {};
  GENOMES.forEach(function (g) { genomeColor[g.genome] = subgroupColor[g[SUBGROUP]]; });

  function hexToRgba(hex, alpha) {
    var r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
    return "rgba(" + r + "," + g + "," + b + "," + alpha + ")";
  }

  // ---- table-aware helpers ---------------------------------------------
  // "table" is "protein" or "cds" everywhere in this file -- which master
  // table a property/sample/stat came from. Kept explicit and threaded
  // through every call rather than inferred, so a property name that
  // happens to exist in both tables (e.g. `length`) is never ambiguous.
  function propsForTable(table) { return table === "cds" ? CDS_PROPERTIES : PROPERTIES; }
  function statsObjFor(table) { return table === "cds" ? DATA.cds_property_stats : DATA.property_stats; }
  function samplesObjFor(table) { return table === "cds" ? DATA.cds_samples : DATA.samples; }

  // ---- grouping-mode helpers ------------------------------------------
  // A "mode" is one of: "genome", PRIMARY, or SUBGROUP -- the three axes a
  // property can be viewed by anywhere in this dashboard.
  function modeOptions() {
    var opts = [{ value: "genome", label: "species (genome)" }];
    opts.push({ value: PRIMARY, label: PRIMARY });
    opts.push({ value: SUBGROUP, label: SUBGROUP });
    return opts;
  }

  function groupsForMode(mode) {
    if (mode === "genome") return GENOME_NAMES;
    if (mode === PRIMARY) return PRIMARY_VALUES;
    if (mode === SUBGROUP) return SUBGROUP_VALUES;
    return [];
  }

  function colorForMode(mode, value) {
    if (mode === "genome") return genomeColor[value];
    if (mode === PRIMARY) return primaryColor[value];
    return subgroupColor[value];
  }

  function statsFor(table, prop, mode, value) {
    var ps = statsObjFor(table)[prop];
    if (!ps) return null;
    if (mode === "genome") return ps.by_genome[value];
    return ps.by_grouping[mode][value];
  }

  var sampleIndexKey = {};
  sampleIndexKey["genome"] = "genome_index";
  sampleIndexKey[PRIMARY] = PRIMARY + "_index";
  sampleIndexKey[SUBGROUP] = SUBGROUP + "_index";

  function sampledValuesForGroup(table, prop, mode, value) {
    var samplesObj = samplesObjFor(table);
    var idxArr = samplesObj[sampleIndexKey[mode]];
    var targetIdx = groupsForMode(mode).indexOf(value);
    var propArr = samplesObj.properties[prop];
    var out = [];
    for (var i = 0; i < idxArr.length; i++) {
      if (idxArr[i] === targetIdx) {
        var v = propArr[i];
        if (v !== null && v !== undefined) out.push(v);
      }
    }
    return out;
  }

  // ---- Gaussian KDE (Plotly has no native KDE trace) -------------------
  function gaussianKDE(values, gridSize) {
    gridSize = gridSize || 200;
    var n = values.length;
    if (n === 0) return { x: [], y: [] };
    var mean = values.reduce(function (a, b) { return a + b; }, 0) / n;
    var variance = n > 1 ? values.reduce(function (a, b) { return a + (b - mean) * (b - mean); }, 0) / (n - 1) : 1;
    var sd = Math.sqrt(variance) || 1;
    var bw = 1.06 * sd * Math.pow(n, -1 / 5) || 1;
    var min = Math.min.apply(null, values), max = Math.max.apply(null, values);
    var pad = (max - min) * 0.15 || 1;
    var x0 = min - pad, x1 = max + pad;
    var xs = new Array(gridSize);
    var ys = new Array(gridSize);
    var norm = 1 / (n * bw * Math.sqrt(2 * Math.PI));
    for (var i = 0; i < gridSize; i++) {
      var x = x0 + (i * (x1 - x0)) / (gridSize - 1);
      var s = 0;
      for (var j = 0; j < n; j++) {
        var u = (x - values[j]) / bw;
        s += Math.exp(-0.5 * u * u);
      }
      xs[i] = x;
      ys[i] = s * norm;
    }
    return { x: xs, y: ys };
  }

  // ---- export helpers ----------------------------------------------------
  function downloadBlob(filename, content, mime) {
    var blob = new Blob([content], { type: mime || "text/csv" });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  function toCSV(columns, rows) {
    var esc = function (v) {
      if (v === null || v === undefined) return "";
      var s = String(v);
      return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    var lines = [columns.map(esc).join(",")];
    rows.forEach(function (row) { lines.push(row.map(esc).join(",")); });
    return lines.join("\n");
  }

  function downloadCSV(filename, columns, rows) {
    downloadBlob(filename, toCSV(columns, rows), "text/csv");
  }

  function downloadPlotPNG(containerId, filename) {
    Plotly.downloadImage(containerId, { format: "png", filename: filename, width: 1100, height: 700 });
  }

  function wireExport(pngBtnId, csvBtnId, plotId, plotFilename, csvFilename, csvFn) {
    var pngBtn = document.getElementById(pngBtnId);
    if (pngBtn) pngBtn.addEventListener("click", function () { downloadPlotPNG(plotId, plotFilename); });
    var csvBtn = document.getElementById(csvBtnId);
    if (csvBtn) csvBtn.addEventListener("click", function () {
      var out = csvFn();
      downloadCSV(csvFilename, out.columns, out.rows);
    });
  }

  // ---- shared render for box / violin / histogram / density ----------
  function renderDistribution(containerId, captionId, table, prop, mode, plotType) {
    var groups = groupsForMode(mode);
    var traces = [];

    if (plotType === "box") {
      groups.forEach(function (g) {
        var s = statsFor(table, prop, mode, g);
        if (!s || !s.n) return;
        traces.push({
          type: "box", name: String(g), x0: String(g),
          q1: [s.q1], median: [s.median], q3: [s.q3],
          lowerfence: [s.min], upperfence: [s.max], mean: [s.mean],
          boxpoints: false, marker: { color: colorForMode(mode, g) }, showlegend: false,
        });
      });
    } else if (plotType === "violin") {
      groups.forEach(function (g) {
        var vals = sampledValuesForGroup(table, prop, mode, g);
        if (!vals.length) return;
        traces.push({
          type: "violin", y: vals, name: String(g), x0: String(g),
          box: { visible: true }, meanline: { visible: true }, points: false,
          marker: { color: colorForMode(mode, g) }, line: { color: colorForMode(mode, g) }, showlegend: false,
        });
      });
    } else if (plotType === "histogram") {
      groups.forEach(function (g) {
        var vals = sampledValuesForGroup(table, prop, mode, g);
        if (!vals.length) return;
        traces.push({
          type: "histogram", x: vals, name: String(g), opacity: 0.55,
          histnorm: "probability density", marker: { color: colorForMode(mode, g) },
        });
      });
    } else if (plotType === "density") {
      groups.forEach(function (g) {
        var vals = sampledValuesForGroup(table, prop, mode, g);
        if (!vals.length) return;
        var kde = gaussianKDE(vals);
        var color = colorForMode(mode, g);
        traces.push({
          type: "scatter", mode: "lines", x: kde.x, y: kde.y, name: String(g),
          line: { color: color, width: 2 }, fill: "tozeroy", fillcolor: hexToRgba(color, 0.15),
        });
      });
    }

    var isDistributional = plotType === "histogram" || plotType === "density";
    var layout = {
      title: { text: prop + " (" + table + ") by " + mode },
      xaxis: { title: { text: isDistributional ? prop : mode } },
      yaxis: { title: { text: isDistributional ? "density" : prop } },
      barmode: plotType === "histogram" ? "overlay" : undefined,
      margin: { t: 46, b: 70, l: 60, r: 20 },
      legend: { orientation: "h", y: -0.25 },
      height: 460,
    };
    Plotly.react(containerId, traces, layout, { responsive: true, displaylogo: false });

    if (captionId) {
      var caption = plotType === "box"
        ? "Exact quartiles computed from the full dataset."
        : "Distribution shown from a " + DATA.meta.sample_per_genome.toLocaleString() + "-record sample per genome (seed=" + DATA.meta.sample_seed + "); summary statistics elsewhere in this dashboard are computed on the full dataset.";
      document.getElementById(captionId).textContent = caption;
    }
  }

  // ---- Property Explorer ----------------------------------------------
  function populateSelect(selectEl, items, valueFn, labelFn, selectValue) {
    selectEl.innerHTML = "";
    items.forEach(function (item) {
      var opt = document.createElement("option");
      opt.value = valueFn(item);
      opt.textContent = labelFn(item);
      selectEl.appendChild(opt);
    });
    if (selectValue !== undefined && items.some(function (i) { return valueFn(i) === selectValue; })) {
      selectEl.value = selectValue;
    }
  }

  var expTable = document.getElementById("expTable");
  var expProperty = document.getElementById("expProperty");
  var expGrouping = document.getElementById("expGrouping");
  var expPlotType = document.getElementById("expPlotType");

  function refreshExpProperties(selectValue) {
    populateSelect(expProperty, propsForTable(expTable.value), function (p) { return p; }, function (p) { return p; }, selectValue);
  }
  refreshExpProperties(DATA.meta.default_property);
  populateSelect(expGrouping, modeOptions(), function (m) { return m.value; }, function (m) { return m.label; });

  function renderExplorer() {
    renderDistribution("expPlot", "expCaption", expTable.value, expProperty.value, expGrouping.value, expPlotType.value);
  }
  expTable.addEventListener("change", function () { refreshExpProperties(); renderExplorer(); });
  [expProperty, expGrouping, expPlotType].forEach(function (el) { el.addEventListener("change", renderExplorer); });

  function jumpToProperty(prop, table) {
    switchSection("explorer");
    expTable.value = table;
    refreshExpProperties(prop);
    renderExplorer();
  }

  wireExport("expExportPng", "expExportCsv", "expPlot", "property_explorer", "property_explorer.csv", function () {
    var table = expTable.value, prop = expProperty.value, mode = expGrouping.value, pt = expPlotType.value;
    var groups = groupsForMode(mode);
    if (pt === "box") {
      var cols = ["group", "n", "min", "q1", "median", "q3", "max", "mean", "std"];
      var rows = groups.map(function (g) {
        var s = statsFor(table, prop, mode, g) || {};
        return [g, s.n, s.min, s.q1, s.median, s.q3, s.max, s.mean, s.std];
      });
      return { columns: cols, rows: rows };
    }
    var rows2 = [];
    groups.forEach(function (g) {
      sampledValuesForGroup(table, prop, mode, g).forEach(function (v) { rows2.push([g, v]); });
    });
    return { columns: ["group", prop], rows: rows2 };
  });

  // ---- Species View -----------------------------------------------------
  var speciesSort = { key: "genome", dir: "asc" };
  var speciesColumns = [];
  function renderSpeciesTable() {
    var table = document.getElementById("speciesTable");
    speciesColumns = ["genome", PRIMARY, SUBGROUP, "n_proteins", "n_cds"].concat(
      PROPERTIES.map(function (p) { return p + " (protein)"; }),
      CDS_PROPERTIES.map(function (p) { return p + " (cds)"; })
    );

    var rows = GENOMES.map(function (g) {
      var row = {
        genome: g.genome,
        n_proteins: g.n_proteins,
        n_cds: (DATA.qc[g.genome] && DATA.qc[g.genome].cds.n_records) || null,
      };
      row[PRIMARY] = g[PRIMARY];
      row[SUBGROUP] = g[SUBGROUP];
      PROPERTIES.forEach(function (p) {
        var s = statsFor("protein", p, "genome", g.genome);
        row[p + " (protein)"] = s ? s.median : null;
      });
      CDS_PROPERTIES.forEach(function (p) {
        var s = statsFor("cds", p, "genome", g.genome);
        row[p + " (cds)"] = s ? s.median : null;
      });
      return row;
    });

    rows.sort(function (a, b) {
      var av = a[speciesSort.key], bv = b[speciesSort.key];
      var cmp;
      if (typeof av === "number" && typeof bv === "number") cmp = av - bv;
      else cmp = String(av).localeCompare(String(bv));
      return speciesSort.dir === "asc" ? cmp : -cmp;
    });

    var thead = "<thead><tr>" + speciesColumns.map(function (c) {
      var cls = c === speciesSort.key ? "sorted " + speciesSort.dir : "";
      return '<th class="' + cls + '" data-key="' + c + '">' + c + "</th>";
    }).join("") + "</tr></thead>";

    var tbody = "<tbody>" + rows.map(function (row) {
      return "<tr>" + speciesColumns.map(function (c) {
        var v = row[c];
        var display = v === null || v === undefined ? "-" : (typeof v === "number" ? (Number.isInteger(v) ? v.toLocaleString() : v.toPrecision(4)) : v);
        return "<td>" + display + "</td>";
      }).join("") + "</tr>";
    }).join("") + "</tbody>";

    table.innerHTML = thead + tbody;
    table._rows = rows;
    table.querySelectorAll("th").forEach(function (th) {
      th.addEventListener("click", function () {
        var key = th.getAttribute("data-key");
        if (speciesSort.key === key) speciesSort.dir = speciesSort.dir === "asc" ? "desc" : "asc";
        else { speciesSort.key = key; speciesSort.dir = "asc"; }
        renderSpeciesTable();
      });
    });
  }

  document.getElementById("spExportCsv").addEventListener("click", function () {
    var table = document.getElementById("speciesTable");
    var rows = (table._rows || []).map(function (row) { return speciesColumns.map(function (c) { return row[c]; }); });
    downloadCSV("species_summary.csv", speciesColumns, rows);
  });

  // ---- Effect Sizes -----------------------------------------------------
  var effGrouping = document.getElementById("effGrouping");

  function effectSizeGroupingOptions() {
    var opts = [{ key: "primary", label: PRIMARY + ": " + PRIMARY_VALUES.join(" vs ") }];
    var seen = {};
    (DATA.effect_sizes[SUBGROUP] || []).forEach(function (r) {
      var k = r.group_a + "|" + r.group_b;
      if (!seen[k]) {
        seen[k] = true;
        opts.push({ key: "sub:" + k, label: SUBGROUP + ": " + r.group_a + " vs " + r.group_b });
      }
    });
    return opts;
  }
  populateSelect(effGrouping, effectSizeGroupingOptions(), function (o) { return o.key; }, function (o) { return o.label; });

  function effectSizeRows(key) {
    if (key === "primary") return DATA.effect_sizes[PRIMARY];
    var pair = key.slice(4).split("|");
    return DATA.effect_sizes[SUBGROUP].filter(function (r) { return r.group_a === pair[0] && r.group_b === pair[1]; });
  }

  var effOrderedRows = [];
  function renderEffectSizes() {
    var rows = effectSizeRows(effGrouping.value).slice();
    rows.sort(function (a, b) { return Math.abs(b.rank_biserial) - Math.abs(a.rank_biserial); });
    effOrderedRows = rows.slice().reverse(); // Plotly draws horizontal bars bottom-to-top
    var labels = effOrderedRows.map(function (r) { return r.property + " (" + r.table + ")"; });
    var values = effOrderedRows.map(function (r) { return r.rank_biserial; });
    var colors = values.map(function (v) { return v < 0 ? "#d62728" : "#1f77b4"; });

    var trace = { type: "bar", orientation: "h", x: values, y: labels, marker: { color: colors } };
    var layout = {
      title: { text: "Effect sizes (rank-biserial correlation)" },
      xaxis: { title: { text: "rank-biserial" }, range: [-1, 1], zeroline: true },
      margin: { l: 230, t: 60, r: 20 },
      height: Math.max(420, 24 * labels.length),
    };
    Plotly.react("effPlot", [trace], layout, { responsive: true, displaylogo: false });

    var plotEl = document.getElementById("effPlot");
    plotEl.removeAllListeners && plotEl.removeAllListeners("plotly_click");
    plotEl.on("plotly_click", function (evt) {
      var point = evt.points[0];
      var row = effOrderedRows[point.pointNumber];
      if (row) jumpToProperty(row.property, row.table);
    });
  }
  effGrouping.addEventListener("change", renderEffectSizes);

  wireExport("effExportPng", "effExportCsv", "effPlot", "effect_sizes", "effect_sizes_comparison.csv", function () {
    var cols = ["property", "table", "group_a", "group_b", "n_a", "n_b", "median_a", "median_b", "p_value", "cles", "rank_biserial"];
    var rows = effOrderedRows.map(function (r) { return cols.map(function (c) { return r[c]; }); });
    return { columns: cols, rows: rows };
  });

  // ---- Sensitivity -------------------------------------------------------
  var sensDisplayed = { propKeys: [], excluded: [], z: [] };
  function renderSensitivity() {
    var rows = DATA.sensitivity.rows;
    var excluded = DATA.sensitivity.excluded_subgroups;
    var byProp = {};
    rows.forEach(function (r) {
      var k = r.property + " (" + r.table + ")";
      byProp[k] = byProp[k] || {};
      byProp[k][r.excluded_subgroup] = r.shrinkage;
    });
    var propKeys = Object.keys(byProp).sort(function (a, b) {
      var maxA = Math.max.apply(null, Object.keys(byProp[a]).map(function (k) { return Math.abs(byProp[a][k]); }));
      var maxB = Math.max.apply(null, Object.keys(byProp[b]).map(function (k) { return Math.abs(byProp[b][k]); }));
      return maxB - maxA;
    });

    // Suppress excluded-subgroup columns that are entirely null (undefined
    // shrinkage for every property -- e.g. dropping a subgroup that leaves
    // zero genomes on one side of the primary grouping) rather than
    // rendering a blank stripe. Detected generically from the data, not
    // from any particular subgroup's name.
    var keptExcluded = excluded.filter(function (s) {
      return propKeys.some(function (k) { return byProp[k][s] !== undefined && byProp[k][s] !== null; });
    });
    var droppedCount = excluded.length - keptExcluded.length;

    var z = propKeys.map(function (k) {
      return keptExcluded.map(function (s) { return byProp[k][s] !== undefined ? byProp[k][s] : null; });
    });
    sensDisplayed = { propKeys: propKeys, excluded: keptExcluded, z: z };
    var maxAbs = 0;
    z.forEach(function (row) { row.forEach(function (v) { if (v !== null) maxAbs = Math.max(maxAbs, Math.abs(v)); }); });
    maxAbs = maxAbs || 1;

    var trace = {
      type: "heatmap", x: keptExcluded, y: propKeys, z: z,
      colorscale: "RdBu", reversescale: true, zmin: -maxAbs, zmax: maxAbs, zmid: 0,
      colorbar: { title: { text: "shrinkage" } }, hoverongaps: false,
    };
    var layout = {
      title: { text: "Leave-one-subgroup-out sensitivity" },
      xaxis: { title: { text: "excluded subgroup" }, type: "category", automargin: true },
      yaxis: { automargin: true },
      margin: { l: 240, t: 46, r: 20, b: 80 },
      height: Math.max(420, 24 * propKeys.length),
    };
    Plotly.react("sensPlot", [trace], layout, { responsive: true, displaylogo: false });

    var noteEl = document.querySelector("#section-sensitivity .note");
    if (noteEl && droppedCount > 0) {
      var base = noteEl.getAttribute("data-base-text") || noteEl.textContent;
      noteEl.setAttribute("data-base-text", base);
      noteEl.textContent = base + " (" + droppedCount + " excluded-subgroup column" + (droppedCount === 1 ? "" : "s") +
        " omitted here: no genomes remain on one side of " + PRIMARY + " once excluded, so shrinkage is undefined.)";
    }
  }

  wireExport("sensExportPng", "sensExportCsv", "sensPlot", "sensitivity", "sensitivity.csv", function () {
    var cols = ["property", "excluded_subgroup", "shrinkage"];
    var rows = [];
    sensDisplayed.propKeys.forEach(function (k) {
      sensDisplayed.excluded.forEach(function (s, j) {
        rows.push([k, s, sensDisplayed.z[sensDisplayed.propKeys.indexOf(k)][j]]);
      });
    });
    return { columns: cols, rows: rows };
  });

  // ---- Codon Usage --------------------------------------------------------
  function renderCodonUsage() {
    var codons = DATA.codon_usage.codons;
    var freqs = DATA.codon_usage.frequencies;
    var z = GENOME_NAMES.map(function (g) { return freqs[g] || codons.map(function () { return null; }); });
    var trace = {
      type: "heatmap", x: codons, y: GENOME_NAMES, z: z,
      colorscale: "YlGnBu", hoverongaps: false,
      colorbar: { title: { text: "% of codons" } },
      hovertemplate: "genome=%{y}<br>codon=%{x}<br>freq=%{z}%<extra></extra>",
    };
    var layout = {
      title: { text: "Codon usage (% frequency per genome)" },
      xaxis: { title: { text: "codon" }, tickangle: -90, tickfont: { size: 9 }, automargin: true },
      yaxis: { automargin: true },
      margin: { l: 140, t: 46, r: 20, b: 90 },
      height: Math.max(420, 28 * GENOME_NAMES.length),
    };
    Plotly.react("codonPlot", [trace], layout, { responsive: true, displaylogo: false });
  }

  wireExport("codonExportPng", "codonExportCsv", "codonPlot", "codon_usage", "codon_usage.csv", function () {
    var codons = DATA.codon_usage.codons;
    var freqs = DATA.codon_usage.frequencies;
    var cols = ["genome"].concat(codons);
    var rows = GENOME_NAMES.map(function (g) { return [g].concat(freqs[g] || []); });
    return { columns: cols, rows: rows };
  });

  // ---- PCA -----------------------------------------------------------------
  var pcaResolution = document.getElementById("pcaResolution");
  var pcaPropertySet = document.getElementById("pcaPropertySet");
  var pcaGrouping = document.getElementById("pcaGrouping");
  var pcaX = document.getElementById("pcaX");
  var pcaY = document.getElementById("pcaY");
  populateSelect(pcaGrouping, modeOptions().filter(function (m) { return m.value !== "genome"; }), function (m) { return m.value; }, function (m) { return m.label; });

  var pcaCurrent = null; // last-rendered PCA block, for export
  var pcaCurrentKind = null; // "species" | "sample"

  function pcaBlockFor(resolution, propSet) {
    if (resolution === "species") return DATA.pca.species[propSet];
    if (resolution === "protein") return DATA.pca.protein_sample;
    return DATA.pca.cds_sample;
  }

  function syncPcaControls() {
    var resolution = pcaResolution.value;
    var forcedSet = resolution === "protein" ? "protein" : resolution === "cds" ? "cds" : null;
    if (forcedSet) pcaPropertySet.value = forcedSet;
    Array.prototype.forEach.call(pcaPropertySet.options, function (opt) {
      opt.disabled = !!forcedSet && opt.value !== forcedSet;
    });
    var block = pcaBlockFor(resolution, pcaPropertySet.value);
    var nComponents = block.explained_variance_ratio.length;
    var pcOptions = [];
    for (var i = 0; i < nComponents; i++) pcOptions.push({ value: i, label: "PC" + (i + 1) + " (" + (block.explained_variance_ratio[i] * 100).toFixed(1) + "% var)" });
    populateSelect(pcaX, pcOptions, function (o) { return o.value; }, function (o) { return o.label; }, 0);
    populateSelect(pcaY, pcOptions, function (o) { return o.value; }, function (o) { return o.label; }, Math.min(1, nComponents - 1));
  }

  function renderLoadings(block, pcX, pcY) {
    var div = document.getElementById("pcaLoadings");
    var entries = block.properties.map(function (p) { return { prop: p, x: block.loadings[p][pcX], y: block.loadings[p][pcY] }; });
    entries.sort(function (a, b) { return (Math.abs(b.x) + Math.abs(b.y)) - (Math.abs(a.x) + Math.abs(a.y)); });
    entries = entries.slice(0, 12);
    var maxAbs = Math.max.apply(null, entries.map(function (e) { return Math.max(Math.abs(e.x), Math.abs(e.y)); }).concat([1e-9]));
    div.innerHTML = entries.map(function (e) {
      var wx = Math.abs(e.x) / maxAbs * 50, wy = Math.abs(e.y) / maxAbs * 50;
      return '<div class="loading-row"><span>' + e.prop + "</span>" +
        '<span class="bar-wrap"><span class="bar" style="width:' + wx.toFixed(1) + '%;background:#1f77b4;"></span></span>' +
        '<span class="bar-wrap"><span class="bar" style="width:' + wy.toFixed(1) + '%;background:#ff7f0e;"></span></span>' +
        "</div>";
    }).join("") + '<p class="caption">Bar length = |loading| on the selected axis (blue = X, orange = Y), top 12 by combined magnitude.</p>';
  }

  function renderPCA() {
    syncPcaControls();
    var resolution = pcaResolution.value, propSet = pcaPropertySet.value;
    var block = pcaBlockFor(resolution, propSet);
    var pcX = Number(pcaX.value), pcY = Number(pcaY.value);
    var mode = pcaGrouping.value;
    pcaCurrent = block; pcaCurrentKind = resolution === "species" ? "species" : "sample";

    var traces = [];
    if (resolution === "species") {
      var order = groupsForMode(mode === PRIMARY ? PRIMARY : SUBGROUP);
      order.forEach(function (gv) {
        var genomesInGroup = GENOMES.filter(function (g) { return g[mode] === gv; }).map(function (g) { return g.genome; });
        var xs = [], ys = [], labels = [];
        genomesInGroup.forEach(function (gname) {
          var s = block.scores[gname];
          if (!s) return;
          xs.push(s[pcX]); ys.push(s[pcY]); labels.push(gname);
        });
        traces.push({
          type: "scatter", mode: "markers+text", x: xs, y: ys, text: labels, textposition: "top center",
          name: String(gv), marker: { color: colorForMode(mode, gv), size: 12 }, textfont: { size: 9 },
        });
      });
    } else {
      var samp = block;
      var idxArr = samp[sampleIndexKey[mode]];
      var groups2 = groupsForMode(mode);
      groups2.forEach(function (gv, gi) {
        var xs = [], ys = [];
        for (var i = 0; i < idxArr.length; i++) {
          if (idxArr[i] === gi) { xs.push(samp.scores[i][pcX]); ys.push(samp.scores[i][pcY]); }
        }
        if (!xs.length) return;
        traces.push({
          type: "scattergl", mode: "markers", x: xs, y: ys, name: String(gv),
          marker: { color: colorForMode(mode, gv), size: 5, opacity: 0.5 },
        });
      });
    }

    var layout = {
      title: { text: (resolution === "species" ? "Per-species" : resolution === "protein" ? "Per-protein (sample)" : "Per-CDS (sample)") + " PCA, colored by " + mode },
      xaxis: { title: { text: "PC" + (pcX + 1) + " (" + (block.explained_variance_ratio[pcX] * 100).toFixed(1) + "% var)" } },
      yaxis: { title: { text: "PC" + (pcY + 1) + " (" + (block.explained_variance_ratio[pcY] * 100).toFixed(1) + "% var)" } },
      margin: { t: 46, b: 60, l: 60, r: 20 },
      legend: { orientation: "h", y: -0.2 },
      height: 500,
    };
    Plotly.react("pcaPlot", traces, layout, { responsive: true, displaylogo: false });
    renderLoadings(block, pcX, pcY);

    document.getElementById("pcaCaption").textContent = resolution === "species"
      ? "Exact PCA on per-genome medians (n=" + Object.keys(block.scores).length + " genomes)."
      : "PCA fit (loadings/explained variance) on the full dataset, exact; points shown are the same " +
        DATA.meta.sample_per_genome.toLocaleString() + "-per-genome sample used elsewhere in this dashboard.";
  }

  [pcaResolution, pcaPropertySet].forEach(function (el) { el.addEventListener("change", renderPCA); });
  [pcaGrouping, pcaX, pcaY].forEach(function (el) { el.addEventListener("change", renderPCA); });

  wireExport("pcaExportPng", "pcaExportCsv", "pcaPlot", "pca", "pca_scores.csv", function () {
    if (pcaCurrentKind === "species") {
      var cols = ["genome"].concat(pcaCurrent.properties.map(function (p, i) { return "PC" + (i + 1); }));
      var rows = Object.keys(pcaCurrent.scores).map(function (g) { return [g].concat(pcaCurrent.scores[g]); });
      return { columns: cols, rows: rows };
    }
    var mode = pcaGrouping.value;
    var idxArr = pcaCurrent[sampleIndexKey[mode]];
    var groups = groupsForMode(mode);
    var nComp = pcaCurrent.explained_variance_ratio.length;
    var cols2 = [mode].concat(Array.from({ length: nComp }, function (_, i) { return "PC" + (i + 1); }));
    var rows2 = pcaCurrent.scores.map(function (row, i) { return [groups[idxArr[i]]].concat(row); });
    return { columns: cols2, rows: rows2 };
  });

  // ---- Clustering: genome x property ---------------------------------------
  function renderClusterGenome() {
    var block = DATA.clustering.genome_property;
    var order = block.genome_order;
    var props = block.properties;
    var z = order.map(function (g) { return props.map(function (p, i) { return block.z[g][i]; }); });

    var heatmap = {
      type: "heatmap", x: props, y: order, z: z, colorscale: "RdBu", reversescale: true, zmid: 0,
      colorbar: { title: { text: "z-score" }, x: 1.02 }, xaxis: "x", yaxis: "y",
    };
    var dendro = block.dendrogram;
    var dendroTraces = [];
    var nLeaves = order.length;
    if (dendro.icoord && dendro.icoord.length) {
      var maxD = Math.max.apply(null, dendro.dcoord.map(function (seg) { return Math.max.apply(null, seg); })) || 1;
      dendro.icoord.forEach(function (xseg, i) {
        var yseg = dendro.dcoord[i];
        dendroTraces.push({
          type: "scatter", mode: "lines", xaxis: "x2", yaxis: "y2",
          x: yseg.map(function (v) { return -v / maxD; }),
          y: xseg.map(function (v) { return (v - 5) / 10; }),
          line: { color: "#6b7280", width: 1 }, hoverinfo: "skip", showlegend: false,
        });
      });
    }

    var layout = {
      title: { text: "Genome x property (z-scored medians)" },
      grid: { rows: 1, columns: 2, pattern: "independent" },
      xaxis: { domain: [0.14, 1], title: { text: "property" }, tickangle: -60, tickfont: { size: 9 }, automargin: true },
      yaxis: { domain: [0, 1], automargin: true },
      xaxis2: { domain: [0, 0.12], visible: false, range: [-1.05, 0.05] },
      yaxis2: { domain: [0, 1], visible: false, range: [-0.5, nLeaves - 0.5] },
      margin: { t: 46, b: 140, l: 20, r: 60 },
      height: Math.max(460, 26 * nLeaves + 160),
      showlegend: false,
    };
    Plotly.react("clusterGenomePlot", [heatmap].concat(dendroTraces), layout, { responsive: true, displaylogo: false });
  }

  wireExport("clusterGenomeExportPng", "clusterGenomeExportCsv", "clusterGenomePlot", "genome_property_clustering", "genome_property_zscores.csv", function () {
    var block = DATA.clustering.genome_property;
    var cols = ["genome"].concat(block.properties);
    var rows = block.genome_order.map(function (g) { return [g].concat(block.z[g]); });
    return { columns: cols, rows: rows };
  });

  // ---- Clustering: property correlation ------------------------------------
  var corrTable = document.getElementById("corrTable");
  var corrCurrent = null;
  function renderCorrelation() {
    var block = DATA.clustering.correlation[corrTable.value];
    corrCurrent = block;
    if (!block.matrix.length) {
      Plotly.react("corrPlot", [], { title: { text: "Not enough properties to correlate" }, height: 300 }, { responsive: true });
      return;
    }
    var trace = {
      type: "heatmap", x: block.properties, y: block.properties, z: block.matrix,
      colorscale: "RdBu", reversescale: true, zmin: -1, zmax: 1, zmid: 0,
      colorbar: { title: { text: "Spearman rho" } },
      hovertemplate: "%{y} vs %{x}<br>rho=%{z}<extra></extra>",
    };
    var layout = {
      title: { text: "Property-property correlation (" + corrTable.value + ", Spearman, clustered)" },
      xaxis: { tickangle: -60, tickfont: { size: 9 }, automargin: true },
      yaxis: { tickfont: { size: 9 }, automargin: true },
      margin: { t: 46, b: 120, l: 140, r: 20 },
      height: Math.max(460, 24 * block.properties.length),
    };
    Plotly.react("corrPlot", [trace], layout, { responsive: true, displaylogo: false });
  }
  corrTable.addEventListener("change", renderCorrelation);

  wireExport("corrExportPng", "corrExportCsv", "corrPlot", "property_correlation", "property_correlation.csv", function () {
    var cols = [""].concat(corrCurrent.properties);
    var rows = corrCurrent.properties.map(function (p, i) { return [p].concat(corrCurrent.matrix[i]); });
    return { columns: cols, rows: rows };
  });

  // ---- Cross-property Scatter -----------------------------------------------
  var scatterTable = document.getElementById("scatterTable");
  var scatterX = document.getElementById("scatterX");
  var scatterY = document.getElementById("scatterY");
  var scatterGrouping = document.getElementById("scatterGrouping");
  var scatterMarginals = document.getElementById("scatterMarginals");
  populateSelect(scatterGrouping, modeOptions(), function (m) { return m.value; }, function (m) { return m.label; });

  function refreshScatterProperties() {
    var props = propsForTable(scatterTable.value);
    populateSelect(scatterX, props, function (p) { return p; }, function (p) { return p; }, props[0]);
    populateSelect(scatterY, props, function (p) { return p; }, function (p) { return p; }, props[1] || props[0]);
  }
  refreshScatterProperties();

  function exactCorrelationFor(table, propX, propY) {
    var block = DATA.clustering.correlation[table];
    var ix = block.properties.indexOf(propX), iy = block.properties.indexOf(propY);
    if (ix === -1 || iy === -1) return null;
    return block.matrix[ix][iy];
  }

  var scatterLastData = null;
  function renderScatter() {
    var table = scatterTable.value, propX = scatterX.value, propY = scatterY.value, mode = scatterGrouping.value;
    var samplesObj = samplesObjFor(table);
    var idxArr = samplesObj[sampleIndexKey[mode]];
    var xArr = samplesObj.properties[propX], yArr = samplesObj.properties[propY];
    var groups = groupsForMode(mode);

    var byGroup = {};
    groups.forEach(function (g) { byGroup[g] = { x: [], y: [] }; });
    for (var i = 0; i < idxArr.length; i++) {
      var xv = xArr[i], yv = yArr[i];
      if (xv === null || yv === null || xv === undefined || yv === undefined) continue;
      var gname = groups[idxArr[i]];
      byGroup[gname].x.push(xv); byGroup[gname].y.push(yv);
    }
    scatterLastData = { table: table, propX: propX, propY: propY, byGroup: byGroup, groups: groups };

    var showMarginals = scatterMarginals.checked;
    var mainXDomain = showMarginals ? [0, 0.78] : [0, 1];
    var mainYDomain = showMarginals ? [0, 0.78] : [0, 1];
    var traces = groups.map(function (g) {
      return {
        type: "scattergl", mode: "markers", x: byGroup[g].x, y: byGroup[g].y, name: String(g),
        marker: { color: colorForMode(mode, g), size: 5, opacity: 0.5 },
        xaxis: "x", yaxis: "y",
      };
    });
    var layout = {
      title: { text: propX + " vs " + propY + " (" + table + ")" },
      xaxis: { title: { text: propX }, domain: mainXDomain },
      yaxis: { title: { text: propY }, domain: mainYDomain },
      margin: { t: 46, b: 60, l: 70, r: 20 },
      legend: { orientation: "h", y: -0.2 },
      height: 520,
    };
    if (showMarginals) {
      // Top marginal (X histograms) shares the main plot's x-axis ('x') but
      // gets its own y-axis ('y3') in the strip above; right marginal (Y
      // histograms) shares the main plot's y-axis ('y') but gets its own
      // x-axis ('x4') in the strip to the right -- standard manual
      // marginal-plot technique in vanilla Plotly.js (no built-in
      // marginal support outside Plotly Express).
      groups.forEach(function (g) {
        traces.push({ type: "histogram", x: byGroup[g].x, xaxis: "x", yaxis: "y3", marker: { color: colorForMode(mode, g) }, opacity: 0.5, showlegend: false });
        traces.push({ type: "histogram", y: byGroup[g].y, xaxis: "x4", yaxis: "y", marker: { color: colorForMode(mode, g) }, opacity: 0.5, showlegend: false });
      });
      layout.yaxis3 = { domain: [0.82, 1], showticklabels: false };
      layout.xaxis4 = { domain: [0.82, 1], showticklabels: false };
      layout.barmode = "overlay";
    }
    Plotly.react("scatterPlot", traces, layout, { responsive: true, displaylogo: false });

    var exact = exactCorrelationFor(table, propX, propY);
    var corrText = exact === null
      ? "Correlation not available (need at least 2 distinct properties in the same table)."
      : "Spearman rho = " + exact + " (exact, computed from the full dataset -- same value as the Clustering view's correlation matrix).";
    document.getElementById("scatterCaption").textContent =
      corrText + " Points shown are the " + DATA.meta.sample_per_genome.toLocaleString() + "-per-genome sample.";
  }
  [scatterX, scatterY, scatterGrouping, scatterMarginals].forEach(function (el) { el.addEventListener("change", renderScatter); });
  scatterTable.addEventListener("change", function () { refreshScatterProperties(); renderScatter(); });

  wireExport("scatterExportPng", "scatterExportCsv", "scatterPlot", "cross_property_scatter", "cross_property_scatter.csv", function () {
    var cols = [scatterGrouping.value, scatterLastData.propX, scatterLastData.propY];
    var rows = [];
    scatterLastData.groups.forEach(function (g) {
      var d = scatterLastData.byGroup[g];
      for (var i = 0; i < d.x.length; i++) rows.push([g, d.x[i], d.y[i]]);
    });
    return { columns: cols, rows: rows };
  });

  // ---- Export section (full summary tables) --------------------------------
  document.getElementById("exportSpeciesSummary").addEventListener("click", function () {
    var cols = ["genome", PRIMARY, SUBGROUP, "n_proteins", "n_cds"];
    PROPERTIES.forEach(function (p) { cols.push(p + "_median", p + "_mean", p + "_std"); });
    CDS_PROPERTIES.forEach(function (p) { cols.push(p + "_cds_median", p + "_cds_mean", p + "_cds_std"); });
    var rows = GENOMES.map(function (g) {
      var row = [g.genome, g[PRIMARY], g[SUBGROUP], g.n_proteins, (DATA.qc[g.genome] && DATA.qc[g.genome].cds.n_records) || null];
      PROPERTIES.forEach(function (p) {
        var s = statsFor("protein", p, "genome", g.genome) || {};
        row.push(s.median, s.mean, s.std);
      });
      CDS_PROPERTIES.forEach(function (p) {
        var s = statsFor("cds", p, "genome", g.genome) || {};
        row.push(s.median, s.mean, s.std);
      });
      return row;
    });
    downloadCSV("species_summary_reconstructed.csv", cols, rows);
  });

  document.getElementById("exportEffectSizes").addEventListener("click", function () {
    var cols = ["grouping", "property", "table", "group_a", "group_b", "n_a", "n_b", "median_a", "median_b", "p_value", "cles", "rank_biserial"];
    var rows = [];
    [PRIMARY, SUBGROUP].forEach(function (grouping) {
      (DATA.effect_sizes[grouping] || []).forEach(function (r) {
        rows.push([grouping, r.property, r.table, r.group_a, r.group_b, r.n_a, r.n_b, r.median_a, r.median_b, r.p_value, r.cles, r.rank_biserial]);
      });
    });
    downloadCSV("effect_sizes_all.csv", cols, rows);
  });

  document.getElementById("exportSensitivity").addEventListener("click", function () {
    var cols = ["excluded_subgroup", "property", "table", "rank_biserial_full", "rank_biserial_excluded", "shrinkage"];
    var rows = DATA.sensitivity.rows.map(function (r) {
      return [r.excluded_subgroup, r.property, r.table, r.rank_biserial_full, r.rank_biserial_excluded, r.shrinkage];
    });
    downloadCSV("sensitivity_leave_one_out.csv", cols, rows);
  });

  // ---- Overview -----------------------------------------------------------
  function renderOverview() {
    document.getElementById("ovIntro").textContent =
      "Cross-genome protein and CDS property comparison across " + DATA.meta.n_genomes + " genomes, grouped by " +
      PRIMARY + " (" + PRIMARY_VALUES.join(", ") + ") and " + SUBGROUP + " (" + SUBGROUP_VALUES.join(", ") + ").";
    document.getElementById("ovNGenomes").textContent = DATA.meta.n_genomes;
    document.getElementById("ovNProteins").textContent = DATA.meta.n_proteins_total.toLocaleString();
    document.getElementById("ovNCds").textContent = DATA.meta.n_cds_total.toLocaleString();
    document.getElementById("ovNProperties").textContent = (PROPERTIES.length + CDS_PROPERTIES.length);

    var groupsHtml = "";
    [PRIMARY, SUBGROUP].forEach(function (col) {
      var values = DATA.grouping_values[col];
      groupsHtml += "<p><strong>" + col + "</strong>: ";
      groupsHtml += values.map(function (v) {
        var color = col === PRIMARY ? primaryColor[v] : subgroupColor[v];
        var count = GENOMES.filter(function (g) { return g[col] === v; }).length;
        return '<span class="legend-swatch" style="background:' + color + '"></span>' + v + " (" + count + " genomes)";
      }).join(" &nbsp; ");
      groupsHtml += "</p>";
    });
    document.getElementById("ovGroups").innerHTML = groupsHtml;

    var qcCols = ["genome", "protein n_records", "n_with_X", "n_with_internal_stop", "cds n_records", "n_non_multiple_of_3", "n_non_ATGC", "record diff", "status"];
    var qcRows = GENOMES.map(function (g) {
      var q = DATA.qc[g.genome] || {};
      var p = q.protein || {}, c = q.cds || {};
      return [
        g.genome, p.n_records, p.n_with_X, p.n_with_internal_stop,
        c.n_records, c.n_non_multiple_of_3, c.n_non_ATGC, q.record_count_diff,
        q.passed ? '<span class="qc-pass">PASSED</span>' : '<span class="qc-fail">FAILED</span>',
      ];
    });
    var thead = "<thead><tr>" + qcCols.map(function (c) { return "<th>" + c + "</th>"; }).join("") + "</tr></thead>";
    var tbody = "<tbody>" + qcRows.map(function (r) { return "<tr>" + r.map(function (v) { return "<td>" + v + "</td>"; }).join("") + "</tr>"; }).join("") + "</tbody>";
    document.getElementById("ovQcTable").innerHTML = thead + tbody;
  }

  // ---- navigation -----------------------------------------------------------
  var sections = {
    overview: renderOverview,
    explorer: renderExplorer,
    species: renderSpeciesTable,
    effects: renderEffectSizes,
    sensitivity: renderSensitivity,
    codon: renderCodonUsage,
    pca: renderPCA,
    clustering: function () { renderClusterGenome(); renderCorrelation(); },
    scatter: renderScatter,
    export: function () {},
  };

  function switchSection(name) {
    document.querySelectorAll(".navbtn").forEach(function (btn) {
      btn.classList.toggle("active", btn.getAttribute("data-section") === name);
    });
    document.querySelectorAll(".section").forEach(function (sec) {
      sec.classList.toggle("active", sec.id === "section-" + name);
    });
    if (sections[name]) sections[name]();
    window.dispatchEvent(new Event("resize"));
  }

  document.querySelectorAll(".navbtn").forEach(function (btn) {
    btn.addEventListener("click", function () { switchSection(btn.getAttribute("data-section")); });
  });

  // sensible defaults so the dashboard looks good instantly on load
  expPlotType.value = "box";
  switchSection("overview");
})();

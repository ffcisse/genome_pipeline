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

  function statsFor(prop, mode, value) {
    var ps = DATA.property_stats[prop];
    if (!ps) return null;
    if (mode === "genome") return ps.by_genome[value];
    return ps.by_grouping[mode][value];
  }

  var sampleIndexKey = {};
  sampleIndexKey["genome"] = "genome_index";
  sampleIndexKey[PRIMARY] = PRIMARY + "_index";
  sampleIndexKey[SUBGROUP] = SUBGROUP + "_index";

  function sampledValuesForGroup(prop, mode, value) {
    var idxArr = DATA.samples[sampleIndexKey[mode]];
    var targetIdx = groupsForMode(mode).indexOf(value);
    var propArr = DATA.samples.properties[prop];
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

  // ---- shared render for box / violin / histogram / density ----------
  function renderDistribution(containerId, captionId, prop, mode, plotType) {
    var groups = groupsForMode(mode);
    var traces = [];

    if (plotType === "box") {
      groups.forEach(function (g) {
        var s = statsFor(prop, mode, g);
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
        var vals = sampledValuesForGroup(prop, mode, g);
        if (!vals.length) return;
        traces.push({
          type: "violin", y: vals, name: String(g), x0: String(g),
          box: { visible: true }, meanline: { visible: true }, points: false,
          marker: { color: colorForMode(mode, g) }, line: { color: colorForMode(mode, g) }, showlegend: false,
        });
      });
    } else if (plotType === "histogram") {
      groups.forEach(function (g) {
        var vals = sampledValuesForGroup(prop, mode, g);
        if (!vals.length) return;
        traces.push({
          type: "histogram", x: vals, name: String(g), opacity: 0.55,
          histnorm: "probability density", marker: { color: colorForMode(mode, g) },
        });
      });
    } else if (plotType === "density") {
      groups.forEach(function (g) {
        var vals = sampledValuesForGroup(prop, mode, g);
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
      title: { text: prop + " by " + mode },
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
        ? "Exact quartiles computed from the full dataset (n=" + DATA.meta.n_proteins_total.toLocaleString() + " proteins)."
        : "Distribution shown from a " + DATA.meta.sample_per_genome.toLocaleString() + "-protein sample per genome (seed=" + DATA.meta.sample_seed + "); summary statistics elsewhere in this dashboard are computed on the full dataset.";
      document.getElementById(captionId).textContent = caption;
    }
  }

  // ---- Property Explorer ----------------------------------------------
  function populateSelect(selectEl, items, valueFn, labelFn) {
    selectEl.innerHTML = "";
    items.forEach(function (item) {
      var opt = document.createElement("option");
      opt.value = valueFn(item);
      opt.textContent = labelFn(item);
      selectEl.appendChild(opt);
    });
  }

  var expProperty = document.getElementById("expProperty");
  var expGrouping = document.getElementById("expGrouping");
  var expPlotType = document.getElementById("expPlotType");

  populateSelect(expProperty, PROPERTIES, function (p) { return p; }, function (p) { return p; });
  populateSelect(expGrouping, modeOptions(), function (m) { return m.value; }, function (m) { return m.label; });

  function renderExplorer() {
    renderDistribution("expPlot", "expCaption", expProperty.value, expGrouping.value, expPlotType.value);
  }
  [expProperty, expGrouping, expPlotType].forEach(function (el) { el.addEventListener("change", renderExplorer); });

  function jumpToProperty(prop) {
    if (PROPERTIES.indexOf(prop) === -1) {
      showToast('"' + prop + '" is a CDS-level property -- its interactive view is planned for Phase 6b.');
      return;
    }
    switchSection("explorer");
    expProperty.value = prop;
    renderExplorer();
  }

  // ---- Species View -----------------------------------------------------
  var spProperty = document.getElementById("spProperty");
  var spPlotType = document.getElementById("spPlotType");
  populateSelect(spProperty, PROPERTIES, function (p) { return p; }, function (p) { return p; });

  function renderSpeciesPlot() {
    renderDistribution("spPlot", "spCaption", spProperty.value, "genome", spPlotType.value);
  }
  [spProperty, spPlotType].forEach(function (el) { el.addEventListener("change", renderSpeciesPlot); });

  var speciesSort = { key: "genome", dir: "asc" };
  function renderSpeciesTable() {
    var table = document.getElementById("speciesTable");
    var columns = ["genome", PRIMARY, SUBGROUP, "n_proteins", "n_cds"].concat(PROPERTIES);

    var rows = GENOMES.map(function (g) {
      var row = {
        genome: g.genome,
        n_proteins: g.n_proteins,
        n_cds: (DATA.qc[g.genome] && DATA.qc[g.genome].cds.n_records) || null,
      };
      row[PRIMARY] = g[PRIMARY];
      row[SUBGROUP] = g[SUBGROUP];
      PROPERTIES.forEach(function (p) {
        var s = statsFor(p, "genome", g.genome);
        row[p] = s ? s.median : null;
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

    var thead = "<thead><tr>" + columns.map(function (c) {
      var cls = c === speciesSort.key ? "sorted " + speciesSort.dir : "";
      return '<th class="' + cls + '" data-key="' + c + '">' + c + "</th>";
    }).join("") + "</tr></thead>";

    var tbody = "<tbody>" + rows.map(function (row) {
      return "<tr>" + columns.map(function (c) {
        var v = row[c];
        var display = v === null || v === undefined ? "-" : (typeof v === "number" ? (Number.isInteger(v) ? v.toLocaleString() : v.toPrecision(4)) : v);
        return "<td>" + display + "</td>";
      }).join("") + "</tr>";
    }).join("") + "</tbody>";

    table.innerHTML = thead + tbody;
    table.querySelectorAll("th").forEach(function (th) {
      th.addEventListener("click", function () {
        var key = th.getAttribute("data-key");
        if (speciesSort.key === key) speciesSort.dir = speciesSort.dir === "asc" ? "desc" : "asc";
        else { speciesSort.key = key; speciesSort.dir = "asc"; }
        renderSpeciesTable();
      });
    });
  }

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
      margin: { l: 230, t: 46, r: 20 },
      height: Math.max(420, 24 * labels.length),
    };
    Plotly.react("effPlot", [trace], layout, { responsive: true, displaylogo: false });

    var plotEl = document.getElementById("effPlot");
    plotEl.removeAllListeners && plotEl.removeAllListeners("plotly_click");
    plotEl.on("plotly_click", function (evt) {
      var point = evt.points[0];
      var row = effOrderedRows[point.pointNumber];
      if (row) jumpToProperty(row.property);
    });
  }
  effGrouping.addEventListener("change", renderEffectSizes);

  // ---- Sensitivity -------------------------------------------------------
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
    var z = propKeys.map(function (k) {
      return excluded.map(function (s) { return byProp[k][s] !== undefined ? byProp[k][s] : null; });
    });
    var maxAbs = 0;
    z.forEach(function (row) { row.forEach(function (v) { if (v !== null) maxAbs = Math.max(maxAbs, Math.abs(v)); }); });
    maxAbs = maxAbs || 1;

    var trace = {
      type: "heatmap", x: excluded, y: propKeys, z: z,
      colorscale: "RdBu", reversescale: true, zmin: -maxAbs, zmax: maxAbs, zmid: 0,
      colorbar: { title: { text: "shrinkage" } },
    };
    var layout = {
      title: { text: "Leave-one-subgroup-out sensitivity" },
      xaxis: { title: { text: "excluded subgroup" } },
      yaxis: { automargin: true },
      margin: { l: 240, t: 46, r: 20 },
      height: Math.max(420, 24 * propKeys.length),
    };
    Plotly.react("sensPlot", [trace], layout, { responsive: true, displaylogo: false });
  }

  // ---- Overview -----------------------------------------------------------
  function renderOverview() {
    document.getElementById("ovIntro").textContent =
      "Cross-genome protein property comparison across " + DATA.meta.n_genomes + " genomes, grouped by " +
      PRIMARY + " (" + PRIMARY_VALUES.join(", ") + ") and " + SUBGROUP + " (" + SUBGROUP_VALUES.join(", ") + "). " +
      DATA.meta.phase6b_note;
    document.getElementById("ovNGenomes").textContent = DATA.meta.n_genomes;
    document.getElementById("ovNProteins").textContent = DATA.meta.n_proteins_total.toLocaleString();
    document.getElementById("ovNCds").textContent = DATA.meta.n_cds_total.toLocaleString();
    document.getElementById("ovNProperties").textContent = PROPERTIES.length;

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

  // ---- toast ---------------------------------------------------------------
  var toastTimer = null;
  function showToast(msg) {
    var el = document.getElementById("toast");
    el.textContent = msg;
    el.classList.add("show");
    if (toastTimer) clearTimeout(toastTimer);
    toastTimer = setTimeout(function () { el.classList.remove("show"); }, 3200);
  }

  // ---- navigation -----------------------------------------------------------
  var sections = {
    overview: renderOverview,
    explorer: renderExplorer,
    species: function () { renderSpeciesTable(); renderSpeciesPlot(); },
    effects: renderEffectSizes,
    sensitivity: renderSensitivity,
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
  spPlotType.value = "box";
  switchSection("overview");
})();

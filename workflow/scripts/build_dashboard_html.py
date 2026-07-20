#!/usr/bin/env python
"""Stage 7 (Phase 6a) -- Assemble the standalone dashboard HTML.

Pure string assembly: inlines the vendored Plotly bundle, the
dashboard_data JSON payload, and this pipeline's custom app JS into
workflow/resources/dashboard/dashboard_template.html's three placeholders
(__PLOTLY_JS__, __DASHBOARD_DATA__, __APP_JS__), and writes the result to
--output. No templating engine, no build step -- the output is one
self-contained .html file with everything already inline, so it opens
directly in a browser (file://) with no server, no install, and no network
requests at all (nothing in it is fetched; it's all embedded at build
time here, not loaded client-side).
"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-json", required=True)
    parser.add_argument("--plotly-js", required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument("--app-js", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.template) as fh:
        template = fh.read()
    with open(args.plotly_js) as fh:
        plotly_js = fh.read()
    with open(args.data_json) as fh:
        data_json = fh.read()
    with open(args.app_js) as fh:
        app_js = fh.read()

    html = template.replace("__PLOTLY_JS__", plotly_js).replace("__DASHBOARD_DATA__", data_json).replace(
        "__APP_JS__", app_js
    )

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(html)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"wrote dashboard -> {args.output} ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()

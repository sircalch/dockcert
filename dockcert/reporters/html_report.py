"""
Interactive HTML report generator for DockCert.
"""

import os
import jinja2
from dockcert.core.scoring import DockingValidationReport

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DockCert Molecular Docking Validation Report</title>
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --pass-color: #10b981;
            --pass-bg: rgba(16, 185, 129, 0.15);
            --warn-color: #f59e0b;
            --warn-bg: rgba(245, 158, 11, 0.15);
            --fail-color: #ef4444;
            --fail-bg: rgba(239, 68, 68, 0.15);
            --accent-blue: #38bdf8;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem 1rem;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            flex-wrap: wrap;
            gap: 1rem;
        }
        .title-group h1 {
            font-size: 2rem;
            font-weight: 700;
            letter-spacing: -0.025em;
        }
        .title-group p {
            color: var(--text-secondary);
            font-size: 0.95rem;
        }
        .status-badge {
            display: inline-flex;
            align-items: center;
            padding: 0.5rem 1.25rem;
            border-radius: 9999px;
            font-weight: 700;
            font-size: 1.1rem;
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }
        .badge-pass { background-color: var(--pass-bg); color: var(--pass-color); border: 1px solid var(--pass-color); }
        .badge-warning { background-color: var(--warn-bg); color: var(--warn-color); border: 1px solid var(--warn-color); }
        .badge-fail { background-color: var(--fail-bg); color: var(--fail-color); border: 1px solid var(--fail-color); }
        
        .grid-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.25rem;
            margin-bottom: 2rem;
        }
        .card {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.25rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
        }
        .card-label {
            font-size: 0.8rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-secondary);
            margin-bottom: 0.4rem;
        }
        .card-value {
            font-size: 1.7rem;
            font-weight: 700;
            color: var(--text-primary);
        }
        .card-subtext {
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 0.25rem;
        }

        .section-title {
            font-size: 1.3rem;
            font-weight: 600;
            margin-bottom: 1rem;
            color: var(--accent-blue);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 2rem;
            font-size: 0.95rem;
            background-color: var(--card-bg);
            border-radius: 0.75rem;
            overflow: hidden;
            border: 1px solid var(--card-border);
        }
        th, td {
            padding: 0.85rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--card-border);
        }
        th {
            background-color: rgba(255, 255, 255, 0.03);
            color: var(--text-secondary);
            font-weight: 600;
            font-size: 0.85rem;
            text-transform: uppercase;
        }
        tr:hover td { background-color: rgba(255, 255, 255, 0.02); }

        .tag {
            display: inline-block;
            padding: 0.2rem 0.6rem;
            border-radius: 0.375rem;
            font-size: 0.75rem;
            font-weight: 700;
        }
        .tag-pass { background-color: var(--pass-bg); color: var(--pass-color); }
        .tag-warning { background-color: var(--warn-bg); color: var(--warn-color); }
        .tag-fail { background-color: var(--fail-bg); color: var(--fail-color); }

        .box {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.25rem;
            margin-bottom: 2rem;
        }
        pre {
            background-color: rgba(0, 0, 0, 0.4);
            padding: 1rem;
            border-radius: 0.5rem;
            overflow-x: auto;
            color: #38bdf8;
            font-family: monospace;
            font-size: 0.9rem;
            border: 1px solid rgba(255, 255, 255, 0.05);
            white-space: pre-wrap;
        }
        .btn-copy {
            background-color: #2563eb;
            color: white;
            border: none;
            padding: 0.4rem 0.8rem;
            border-radius: 0.375rem;
            cursor: pointer;
            font-size: 0.8rem;
            margin-top: 0.5rem;
        }
        .btn-copy:hover { background-color: #1d4ed8; }

        footer {
            text-align: center;
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 3rem;
            padding-top: 1.5rem;
            border-top: 1px solid var(--card-border);
        }
    </style>
</head>
<body>
    <div class="container">
        <header class="header">
            <div class="title-group">
                <h1>DockCert Validation Dashboard</h1>
                <p>Automated Quality Assurance for Molecular Docking and Virtual Screening</p>
            </div>
            <div>
                <span class="status-badge badge-{{ report.overall_status.lower() }}">
                    {{ report.overall_status }}
                </span>
            </div>
        </header>

        <div class="grid-cards">
            <div class="card">
                <div class="card-label">Validation Status</div>
                <div class="card-value" style="font-size: 1.2rem;">{{ report.validation_score }}</div>
                <div class="card-subtext">{{ report.dataset_summary.n_total if report.dataset_summary else 0 }} compounds screened</div>
            </div>
            {% if report.enrichment_metrics.roc_auc %}
            <div class="card">
                <div class="card-label">ROC-AUC (95% CI)</div>
                <div class="card-value">{{ "%.3f"|format(report.enrichment_metrics.roc_auc.value) }}</div>
                <div class="card-subtext">
                    {% if report.enrichment_metrics.roc_auc.ci_lower_95 %}
                    [{{ "%.2f"|format(report.enrichment_metrics.roc_auc.ci_lower_95) }}, {{ "%.2f"|format(report.enrichment_metrics.roc_auc.ci_upper_95) }}]
                    {% else %}
                    Stratified Bootstrap
                    {% endif %}
                </div>
            </div>
            {% endif %}
            {% if report.enrichment_metrics.bedroc_20 %}
            <div class="card">
                <div class="card-label">BEDROC (&alpha;=20.0)</div>
                <div class="card-value">{{ "%.3f"|format(report.enrichment_metrics.bedroc_20.value) }}</div>
                <div class="card-subtext">Early enrichment (top 8%)</div>
            </div>
            {% endif %}
            {% if report.enrichment_metrics.ef_1pct %}
            <div class="card">
                <div class="card-label">Enrichment Factor (1%)</div>
                <div class="card-value">{{ "%.1f"|format(report.enrichment_metrics.ef_1pct.value) }} &times;</div>
                <div class="card-subtext">Top 1% active recovery</div>
            </div>
            {% endif %}
            {% if report.redocking_result %}
            <div class="card">
                <div class="card-label">Redocking RMSD</div>
                <div class="card-value">{{ "%.2f"|format(report.redocking_result.min_rmsd) }} <span style="font-size: 0.9rem; color: var(--text-secondary)">&Aring;</span></div>
                <div class="card-subtext">Status: {{ report.redocking_result.status }}</div>
            </div>
            {% endif %}
        </div>

        <h2 class="section-title">Enrichment & Validation Metrics</h2>
        <table>
            <thead>
                <tr>
                    <th>Metric</th>
                    <th>Observed Value</th>
                    <th>95% Bootstrap CI</th>
                    <th>Pass Threshold</th>
                    <th>Status</th>
                    <th>Diagnostic Finding</th>
                </tr>
            </thead>
            <tbody>
                {% for key, item in report.enrichment_metrics.items() %}
                <tr>
                    <td><strong>{{ item.name }}</strong></td>
                    <td>{{ "%.3f"|format(item.value) }}</td>
                    <td>
                        {% if item.ci_lower_95 is not none %}
                        [{{ "%.3f"|format(item.ci_lower_95) }}, {{ "%.3f"|format(item.ci_upper_95) }}]
                        {% else %}
                        N/A
                        {% endif %}
                    </td>
                    <td>&ge; {{ "%.2f"|format(item.threshold_pass) }}</td>
                    <td><span class="tag tag-{{ item.status.lower() }}">{{ item.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ item.message }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>

        {% if report.recommendations %}
        <div class="box" style="border-left: 4px solid var(--warn-color);">
            <h3 style="color: var(--warn-color); margin-bottom: 0.5rem;">Diagnostic Notes & Recommendations</h3>
            <ul style="padding-left: 1.25rem;">
                {% for rec in report.recommendations %}
                <li style="margin-bottom: 0.25rem; color: var(--text-secondary);">{{ rec }}</li>
                {% endfor %}
            </ul>
        </div>
        {% endif %}

        <h2 class="section-title">Methods Section Snippet (Publication Ready)</h2>
        <div class="box">
            <pre id="methodsSnippet">{{ methods_text }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('methodsSnippet')">Copy Methods Paragraph</button>
        </div>

        <h2 class="section-title">BibTeX Citation</h2>
        <div class="box">
            <pre id="bibSnippet">{{ citation_bib }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('bibSnippet')">Copy BibTeX</button>
        </div>

        <footer>
            Generated automatically by <strong>DockCert v1.0.0</strong> &bull; Molecular Docking Validation & Statistical Quality Toolkit &bull; Monreal-Hernández, 2026.
        </footer>
    </div>

    <script>
        function copyToClipboard(elementId) {
            const text = document.getElementById(elementId).innerText;
            navigator.clipboard.writeText(text).then(() => {
                alert('Copied to clipboard!');
            }).catch(err => {
                console.error('Error copying: ', err);
            });
        }
    </script>
</body>
</html>
"""


def generate_docking_html_report(
    report: DockingValidationReport,
    output_path: str,
    methods_text: str = "",
    citation_bib: str = ""
) -> str:
    """
    Renders the HTML report template and writes it to disk.
    """
    template = jinja2.Template(HTML_TEMPLATE)
    rendered = template.render(
        report=report,
        methods_text=methods_text,
        citation_bib=citation_bib
    )
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered)
    return output_path

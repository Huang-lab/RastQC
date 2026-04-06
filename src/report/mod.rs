use crate::modules::{QCModule, QCResult};
use crate::FileSummary;
use std::io::Write;
use serde_json;

pub fn generate_html_report(filename: &str, modules: &[Box<dyn QCModule>]) -> String {
    let mut html = String::with_capacity(64 * 1024);

    html.push_str("<!DOCTYPE html>\n<html>\n<head>\n");
    html.push_str(&format!(
        "<title>RastQC Report: {}</title>\n",
        html_escape(filename)
    ));
    html.push_str("<meta charset=\"utf-8\">\n");
    html.push_str("<style>\n");
    html.push_str(CSS_STYLES);
    html.push_str("</style>\n");
    html.push_str("</head>\n<body>\n");

    // Header
    html.push_str("<div id=\"header\">\n");
    html.push_str("<h1>RastQC Report</h1>\n");
    html.push_str(&format!(
        "<h2>{}</h2>\n",
        html_escape(filename)
    ));
    html.push_str("</div>\n");

    // Summary sidebar
    html.push_str("<div id=\"sidebar\">\n");
    html.push_str("<h3>Summary</h3>\n");
    html.push_str("<ul>\n");
    for module in modules {
        let icon = match module.result() {
            QCResult::Pass => "&#x2714;",
            QCResult::Warn => "&#x26A0;",
            QCResult::Fail => "&#x2718;",
            QCResult::NotRun => "&#x2796;",
        };
        let class = module.result().icon();
        let id = module.name().replace(' ', "_").to_lowercase();
        html.push_str(&format!(
            "<li class=\"{}\"><a href=\"#{}\">{} {}</a></li>\n",
            class,
            id,
            icon,
            html_escape(module.name())
        ));
    }
    html.push_str("</ul>\n</div>\n");

    // Main content
    html.push_str("<div id=\"main\">\n");

    for module in modules {
        let id = module.name().replace(' ', "_").to_lowercase();
        let status_class = module.result().icon();

        html.push_str(&format!(
            "<div class=\"module {}\" id=\"{}\">\n",
            status_class, id
        ));
        html.push_str(&format!(
            "<h2>{} <span class=\"status {}\">[{}]</span></h2>\n",
            html_escape(module.name()),
            status_class,
            module.result().label()
        ));

        if module.has_chart() {
            let svg = module.svg_chart();
            if !svg.is_empty() {
                html.push_str("<div class=\"chart\">\n");
                html.push_str(&svg);
                html.push_str("\n</div>\n");
            }
        }

        // For table-based modules (Basic Stats, Overrepresented Seqs, Kmer Content)
        let text = module.text_data();
        if !module.has_chart() {
            html.push_str(&text_data_to_html_table(&text));
        }

        html.push_str("</div>\n");
    }

    html.push_str("</div>\n");

    // Footer
    html.push_str("<div id=\"footer\">\n");
    html.push_str("<p>Produced by <strong>RastQC</strong> v0.1.0</p>\n");
    html.push_str("</div>\n");

    html.push_str("</body>\n</html>");
    html
}

fn text_data_to_html_table(text: &str) -> String {
    let mut html = String::new();
    html.push_str("<table class=\"data-table\">\n");

    let mut in_module = false;
    for line in text.lines() {
        if line.starts_with(">>") {
            if line.starts_with(">>END_MODULE") {
                break;
            }
            in_module = true;
            continue;
        }
        if !in_module {
            continue;
        }

        let tag = if line.starts_with('#') { "th" } else { "td" };
        let line = line.trim_start_matches('#');

        html.push_str("<tr>");
        for cell in line.split('\t') {
            html.push_str(&format!("<{tag}>{}</{tag}>", html_escape(cell)));
        }
        html.push_str("</tr>\n");
    }

    html.push_str("</table>\n");
    html
}

pub fn generate_text_data(filename: &str, modules: &[Box<dyn QCModule>]) -> String {
    let mut text = String::with_capacity(32 * 1024);
    text.push_str("##FastQC\t0.12.0\n");

    for module in modules {
        text.push_str(&module.text_data());
    }

    // Inject filename into Basic Statistics (replaces {} placeholder)
    text = text.replacen("Filename\t{}", &format!("Filename\t{}", filename), 1);

    text
}

/// Generate summary.txt content: one line per module with status, name, filename.
pub fn generate_summary_txt(filename: &str, modules: &[Box<dyn QCModule>]) -> String {
    let mut out = String::new();
    for module in modules {
        out.push_str(&format!(
            "{}\t{}\t{}\n",
            module.result().label(),
            module.name(),
            filename
        ));
    }
    out
}

/// Generate native MultiQC JSON output for a single file's QC modules.
/// This produces a JSON object keyed by module name with structured data,
/// eliminating the need for MultiQC to parse fastqc_data.txt.
pub fn generate_multiqc_json(filename: &str, modules: &[Box<dyn QCModule>]) -> String {
    let mut top = serde_json::Map::new();

    // Report metadata
    top.insert(
        "report_metadata".to_string(),
        serde_json::json!({
            "tool": "RastQC",
            "tool_version": "0.1.0",
            "filename": filename,
            "format_version": "1.0"
        }),
    );

    // Module data
    let mut module_data = serde_json::Map::new();
    for module in modules {
        let key = module.name().replace(' ', "_").to_lowercase();
        module_data.insert(key, module.json_data());
    }
    top.insert("modules".to_string(), serde_json::Value::Object(module_data));

    serde_json::to_string_pretty(&serde_json::Value::Object(top)).unwrap_or_default()
}

pub fn write_zip_archive(
    path: &std::path::Path,
    stem: &str,
    html: &str,
    text: &str,
    summary_txt: &str,
) -> anyhow::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut zip = zip::ZipWriter::new(file);

    let options = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);

    let prefix = format!("{}_fastqc", stem);

    zip.start_file(format!("{}/fastqc_report.html", prefix), options)?;
    zip.write_all(html.as_bytes())?;

    zip.start_file(format!("{}/fastqc_data.txt", prefix), options)?;
    zip.write_all(text.as_bytes())?;

    zip.start_file(format!("{}/summary.txt", prefix), options)?;
    zip.write_all(summary_txt.as_bytes())?;

    zip.finish()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Multi-file summary outputs
// ---------------------------------------------------------------------------

/// Generate a TSV summary: rows = files, columns = modules.
pub fn generate_summary_tsv(summaries: &[FileSummary]) -> String {
    if summaries.is_empty() {
        return String::new();
    }

    let mut tsv = String::with_capacity(4096);

    // Header row
    tsv.push_str("Sample");
    for (name, _) in &summaries[0].module_results {
        tsv.push('\t');
        tsv.push_str(name);
    }
    tsv.push_str("\tTotal Sequences\n");

    // Data rows
    for summary in summaries {
        tsv.push_str(&summary.filename);
        for (_, result) in &summary.module_results {
            tsv.push('\t');
            tsv.push_str(result.label());
        }
        tsv.push('\t');
        tsv.push_str(&summary.total_sequences.to_string());
        tsv.push('\n');
    }

    tsv
}

/// Generate an HTML dashboard summarising all files.
pub fn generate_summary_html(summaries: &[FileSummary]) -> String {
    if summaries.is_empty() {
        return String::from("<html><body><p>No files processed.</p></body></html>");
    }

    let mut html = String::with_capacity(16 * 1024);

    html.push_str("<!DOCTYPE html>\n<html>\n<head>\n");
    html.push_str("<title>RastQC Summary</title>\n");
    html.push_str("<meta charset=\"utf-8\">\n");
    html.push_str("<style>\n");
    html.push_str(SUMMARY_CSS);
    html.push_str("</style>\n");
    html.push_str("</head>\n<body>\n");

    // Header
    html.push_str("<div class=\"header\">\n");
    html.push_str("<h1>RastQC Summary</h1>\n");
    html.push_str(&format!(
        "<p>{} samples analysed</p>\n",
        summaries.len()
    ));
    html.push_str("</div>\n");

    // Tally
    let (total_pass, total_warn, total_fail) = count_results(summaries);
    html.push_str("<div class=\"tally\">\n");
    html.push_str(&format!(
        "<span class=\"tally-pass\">{} pass</span>\n",
        total_pass
    ));
    html.push_str(&format!(
        "<span class=\"tally-warn\">{} warn</span>\n",
        total_warn
    ));
    html.push_str(&format!(
        "<span class=\"tally-fail\">{} fail</span>\n",
        total_fail
    ));
    html.push_str("</div>\n");

    // Matrix table
    html.push_str("<div class=\"table-wrap\">\n");
    html.push_str("<table>\n<thead>\n<tr>\n<th>Sample</th>\n");

    // Column headers = module names (abbreviated)
    for (name, _) in &summaries[0].module_results {
        let abbrev = abbreviate_module(name);
        html.push_str(&format!(
            "<th title=\"{}\">{}</th>\n",
            html_escape(name),
            html_escape(&abbrev)
        ));
    }
    html.push_str("<th>Seqs</th>\n</tr>\n</thead>\n<tbody>\n");

    // One row per file
    for summary in summaries {
        // Determine row class from worst result
        let worst = worst_result(&summary.module_results);
        let row_class = worst.icon();

        html.push_str(&format!("<tr class=\"{}\">\n", row_class));
        html.push_str(&format!(
            "<td class=\"sample\"><a href=\"{}\">{}</a></td>\n",
            html_escape(&summary.report_path),
            html_escape(&summary.filename)
        ));

        for (_, result) in &summary.module_results {
            let cls = result.icon();
            let sym = match result {
                QCResult::Pass => "&#x2714;",
                QCResult::Warn => "&#x26A0;",
                QCResult::Fail => "&#x2718;",
                QCResult::NotRun => "-",
            };
            html.push_str(&format!("<td class=\"cell {}\">{}</td>\n", cls, sym));
        }

        html.push_str(&format!(
            "<td class=\"seqs\">{}</td>\n",
            format_count(summary.total_sequences)
        ));
        html.push_str("</tr>\n");
    }

    html.push_str("</tbody>\n</table>\n</div>\n");

    // Footer
    html.push_str("<div class=\"footer\"><p>Produced by <strong>RastQC</strong> v0.1.0</p></div>\n");
    html.push_str("</body>\n</html>");
    html
}

fn count_results(summaries: &[FileSummary]) -> (usize, usize, usize) {
    let mut pass = 0;
    let mut warn = 0;
    let mut fail = 0;
    for s in summaries {
        for (_, r) in &s.module_results {
            match r {
                QCResult::Pass => pass += 1,
                QCResult::Warn => warn += 1,
                QCResult::Fail => fail += 1,
                QCResult::NotRun => {}
            }
        }
    }
    (pass, warn, fail)
}

fn worst_result(results: &[(String, QCResult)]) -> QCResult {
    let mut worst = QCResult::Pass;
    for (_, r) in results {
        match r {
            QCResult::Fail => return QCResult::Fail,
            QCResult::Warn => worst = QCResult::Warn,
            _ => {}
        }
    }
    worst
}

fn abbreviate_module(name: &str) -> String {
    match name {
        "Basic Statistics" => "Basic".into(),
        "Per base sequence quality" => "Base Qual".into(),
        "Per tile sequence quality" => "Tile Qual".into(),
        "Per sequence quality scores" => "Seq Qual".into(),
        "Per base sequence content" => "Base Content".into(),
        "Per sequence GC content" => "GC Content".into(),
        "Per base N content" => "N Content".into(),
        "Sequence Length Distribution" => "Seq Length".into(),
        "Sequence Duplication Levels" => "Duplication".into(),
        "Overrepresented sequences" => "Overrep".into(),
        "Adapter Content" => "Adapters".into(),
        "Kmer Content" => "Kmers".into(),
        "Read Length N50" => "N50".into(),
        "Quality Stratified Length" => "Q-Strat Len".into(),
        "Homopolymer Content" => "Homopolymer".into(),
        other => other.into(),
    }
}

fn format_count(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}K", n as f64 / 1_000.0)
    } else {
        format!("{}", n)
    }
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

const CSS_STYLES: &str = r##"
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif; color: #333; line-height: 1.5; }

#header {
    background: linear-gradient(135deg, #1a5276, #2980b9);
    color: white;
    padding: 20px 30px;
}
#header h1 { font-size: 24px; font-weight: 600; }
#header h2 { font-size: 14px; font-weight: 400; opacity: 0.9; margin-top: 4px; }

#sidebar {
    position: fixed;
    left: 0;
    top: 80px;
    width: 260px;
    height: calc(100vh - 80px);
    overflow-y: auto;
    background: #f8f9fa;
    border-right: 1px solid #dee2e6;
    padding: 15px;
}
#sidebar h3 { font-size: 14px; margin-bottom: 10px; color: #555; text-transform: uppercase; letter-spacing: 1px; }
#sidebar ul { list-style: none; }
#sidebar li { margin: 4px 0; }
#sidebar a { text-decoration: none; color: #333; font-size: 13px; display: block; padding: 4px 8px; border-radius: 4px; }
#sidebar a:hover { background: #e9ecef; }
#sidebar .pass a::before { content: ""; display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: #27ae60; margin-right: 8px; }
#sidebar .warn a::before { content: ""; display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: #f39c12; margin-right: 8px; }
#sidebar .fail a::before { content: ""; display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: #e74c3c; margin-right: 8px; }
#sidebar .na a::before { content: ""; display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: #95a5a6; margin-right: 8px; }

#main {
    margin-left: 280px;
    padding: 20px 30px;
}

.module {
    margin-bottom: 30px;
    padding: 20px;
    border: 1px solid #dee2e6;
    border-radius: 8px;
    background: white;
}
.module h2 { font-size: 18px; margin-bottom: 15px; padding-bottom: 8px; border-bottom: 2px solid #eee; }
.module.pass h2 { border-bottom-color: #27ae60; }
.module.warn h2 { border-bottom-color: #f39c12; }
.module.fail h2 { border-bottom-color: #e74c3c; }

.status { font-size: 12px; font-weight: 600; padding: 2px 8px; border-radius: 3px; }
.status.pass { color: #27ae60; background: #eafaf1; }
.status.warn { color: #f39c12; background: #fef9e7; }
.status.fail { color: #e74c3c; background: #fdedec; }

.chart { margin: 10px 0; max-width: 100%; }
.chart svg { width: 100%; height: auto; max-height: 450px; }

.data-table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }
.data-table th, .data-table td { padding: 6px 12px; border: 1px solid #dee2e6; text-align: left; }
.data-table th { background: #f8f9fa; font-weight: 600; }
.data-table tr:hover { background: #f8f9fa; }

#footer {
    margin-left: 280px;
    padding: 15px 30px;
    border-top: 1px solid #dee2e6;
    color: #888;
    font-size: 12px;
}

@media print {
    #sidebar { display: none; }
    #main, #footer { margin-left: 0; }
}
"##;

const SUMMARY_CSS: &str = r##"
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
    color: #333; line-height: 1.5; background: #f5f6fa;
}

.header {
    background: linear-gradient(135deg, #1a5276, #2980b9);
    color: white; padding: 20px 30px;
}
.header h1 { font-size: 24px; }
.header p { opacity: 0.9; margin-top: 4px; }

.tally {
    padding: 15px 30px; background: white;
    border-bottom: 1px solid #dee2e6;
    display: flex; gap: 20px;
}
.tally span {
    font-size: 14px; font-weight: 600; padding: 4px 12px;
    border-radius: 4px;
}
.tally-pass { color: #27ae60; background: #eafaf1; }
.tally-warn { color: #f39c12; background: #fef9e7; }
.tally-fail { color: #e74c3c; background: #fdedec; }

.table-wrap {
    padding: 20px 30px; overflow-x: auto;
}

table {
    border-collapse: collapse; width: 100%;
    background: white; border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    font-size: 13px;
}
thead th {
    background: #f8f9fa; padding: 8px 6px;
    border: 1px solid #dee2e6; font-weight: 600;
    text-align: center; white-space: nowrap;
    position: sticky; top: 0; z-index: 1;
}
thead th:first-child { text-align: left; min-width: 200px; }

td { padding: 6px; border: 1px solid #dee2e6; text-align: center; }
td.sample { text-align: left; font-weight: 500; white-space: nowrap; }
td.sample a { color: #2980b9; text-decoration: none; }
td.sample a:hover { text-decoration: underline; }
td.seqs { color: #888; font-size: 12px; white-space: nowrap; }

.cell.pass { background: #eafaf1; color: #27ae60; }
.cell.warn { background: #fef9e7; color: #f39c12; }
.cell.fail { background: #fdedec; color: #e74c3c; font-weight: bold; }
.cell.na   { background: #f8f9fa; color: #95a5a6; }

tr:hover td { background-color: #f0f4ff; }
tr:hover td.cell.pass { background: #d5f5e3; }
tr:hover td.cell.warn { background: #fdebd0; }
tr:hover td.cell.fail { background: #fadbd8; }

.footer {
    padding: 15px 30px; color: #888; font-size: 12px;
    border-top: 1px solid #dee2e6;
}
"##;

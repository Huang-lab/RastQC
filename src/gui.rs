use anyhow::Result;
use std::path::Path;
use tiny_http::{Header, Response, Server, StatusCode};

/// Start a local web server that serves RastQC reports and provides
/// a GUI for viewing results and selecting files for analysis.
pub fn start_server(outdir: &Path, port: u16) -> Result<()> {
    let addr = format!("0.0.0.0:{}", port);
    let server = Server::http(&addr)
        .map_err(|e| anyhow::anyhow!("Failed to start server on {}: {}", addr, e))?;

    eprintln!("RastQC GUI server running at http://localhost:{}", port);
    eprintln!("Press Ctrl+C to stop");

    // Try to open in browser
    let url = format!("http://localhost:{}", port);
    if open::that(&url).is_err() {
        eprintln!("Open {} in your browser", url);
    }

    let outdir = outdir.to_path_buf();

    for request in server.incoming_requests() {
        let url_path = request.url().to_string();
        let url_path = url_path.split('?').next().unwrap_or("/");

        let response = match url_path {
            "/" => serve_index(&outdir),
            _ => {
                // Serve static files from outdir
                let file_path = url_path.trim_start_matches('/');
                // Prevent directory traversal
                if file_path.contains("..") {
                    Response::from_string("Forbidden")
                        .with_status_code(StatusCode(403))
                        .with_header(content_type("text/plain"))
                } else {
                    serve_file(&outdir, file_path)
                }
            }
        };

        let _ = request.respond(response);
    }

    Ok(())
}

fn serve_index(outdir: &Path) -> Response<std::io::Cursor<Vec<u8>>> {
    let mut reports: Vec<(String, String)> = Vec::new();

    // Find HTML reports
    if let Ok(entries) = std::fs::read_dir(outdir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            if name.ends_with("_fastqc.html") {
                let display = name.trim_end_matches("_fastqc.html").to_string();
                reports.push((name, display));
            }
        }
    }
    reports.sort_by(|a, b| a.1.cmp(&b.1));

    // Check for summary
    let has_summary = outdir.join("summary.html").exists();

    let mut html = String::with_capacity(8192);
    html.push_str("<!DOCTYPE html>\n<html>\n<head>\n");
    html.push_str("<title>RastQC - Quality Control Reports</title>\n");
    html.push_str("<meta charset=\"utf-8\">\n");
    html.push_str("<style>\n");
    html.push_str(INDEX_CSS);
    html.push_str("</style>\n");
    html.push_str("</head>\n<body>\n");

    html.push_str("<div class=\"header\">\n");
    html.push_str("<h1>RastQC Report Viewer</h1>\n");
    html.push_str(&format!("<p>{} reports available</p>\n", reports.len()));
    html.push_str("</div>\n");

    html.push_str("<div class=\"content\">\n");

    if has_summary {
        html.push_str("<div class=\"summary-link\">\n");
        html.push_str("<a href=\"/summary.html\" class=\"btn\">View Multi-Sample Summary</a>\n");
        html.push_str("</div>\n");
    }

    if reports.is_empty() {
        html.push_str("<p class=\"empty\">No reports found. Run <code>rastqc</code> to generate reports, then refresh this page.</p>\n");
    } else {
        html.push_str("<div class=\"reports\">\n");
        html.push_str("<h2>Individual Reports</h2>\n");
        html.push_str("<ul>\n");
        for (filename, display) in &reports {
            html.push_str(&format!(
                "<li><a href=\"/{}\">{}</a></li>\n",
                html_escape(filename),
                html_escape(display)
            ));
        }
        html.push_str("</ul>\n");
        html.push_str("</div>\n");
    }

    html.push_str("</div>\n");

    html.push_str("<div class=\"footer\">\n");
    html.push_str("<p>Powered by <strong>RastQC</strong> v0.1.0</p>\n");
    html.push_str("</div>\n");
    html.push_str("</body>\n</html>");

    Response::from_data(html.into_bytes()).with_header(content_type("text/html; charset=utf-8"))
}

fn serve_file(outdir: &Path, relative: &str) -> Response<std::io::Cursor<Vec<u8>>> {
    let path = outdir.join(relative);

    match std::fs::read(&path) {
        Ok(data) => {
            let ct = match path.extension().and_then(|e| e.to_str()) {
                Some("html") => "text/html; charset=utf-8",
                Some("txt") | Some("tsv") => "text/plain; charset=utf-8",
                Some("json") => "application/json; charset=utf-8",
                Some("zip") => "application/zip",
                Some("css") => "text/css; charset=utf-8",
                Some("js") => "application/javascript; charset=utf-8",
                Some("svg") => "image/svg+xml",
                Some("png") => "image/png",
                _ => "application/octet-stream",
            };
            Response::from_data(data).with_header(content_type(ct))
        }
        Err(_) => Response::from_string("Not Found")
            .with_status_code(StatusCode(404))
            .with_header(content_type("text/plain")),
    }
}

fn content_type(ct: &str) -> Header {
    Header::from_bytes("Content-Type", ct).unwrap()
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

const INDEX_CSS: &str = r##"
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
    color: #333; line-height: 1.6; background: #f5f6fa;
}
.header {
    background: linear-gradient(135deg, #1a5276, #2980b9);
    color: white; padding: 30px 40px;
}
.header h1 { font-size: 28px; }
.header p { opacity: 0.9; margin-top: 6px; }
.content { max-width: 900px; margin: 30px auto; padding: 0 20px; }
.summary-link { margin-bottom: 20px; }
.btn {
    display: inline-block; padding: 10px 24px;
    background: #2980b9; color: white; text-decoration: none;
    border-radius: 6px; font-weight: 600; font-size: 14px;
}
.btn:hover { background: #1a5276; }
.reports { background: white; border-radius: 8px; padding: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
.reports h2 { font-size: 18px; margin-bottom: 12px; color: #1a5276; }
.reports ul { list-style: none; }
.reports li {
    padding: 10px 12px; border-bottom: 1px solid #eee;
}
.reports li:last-child { border-bottom: none; }
.reports a { color: #2980b9; text-decoration: none; font-size: 14px; }
.reports a:hover { text-decoration: underline; }
.empty { color: #888; font-style: italic; padding: 20px 0; }
.footer { text-align: center; padding: 20px; color: #888; font-size: 12px; }
"##;

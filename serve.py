"""Local development server: like `python3 -m http.server`, but tells the browser to revalidate every
file (Cache-Control: no-cache), so an edited ES module is never mixed with a stale cached one.
Unchanged files still come back as a cheap 304.

    python3 serve.py [port]      # default 8000, then open http://localhost:8000
"""
import http.server
import sys


class NoCacheHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f"serving on http://localhost:{port}")
    http.server.ThreadingHTTPServer(("", port), NoCacheHandler).serve_forever()

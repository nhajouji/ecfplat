"""Inject social-share (Open Graph / Twitter) meta tags into Streamlit.

Streamlit has no first-class way to set link-preview metadata: the HTML
shell it serves is `static/index.html` inside the installed package, and
crawlers never run the app's JavaScript. So we patch that file once at
deploy time (see the startCommand in render.yaml). Idempotent — a marker
comment keeps repeat runs from double-inserting.

The image itself is served by the app via enableStaticServing
(.streamlit/config.toml) from ./static, at /app/static/og_image.jpg.
"""

from pathlib import Path

import streamlit

SITE = "https://elliptic-curves.info"
TITLE = "Elliptic Curves"
DESC = ("Pictures, tools, and theory — from curves over finite fields to "
        "the lattices that draw them. By Nadir Hajouji and Steve Trettel.")
IMAGE = f"{SITE}/app/static/og_image.jpg"
MARK = "<!-- ecfplat og tags -->"

TAGS = f"""{MARK}
<meta property="og:type" content="website">
<meta property="og:site_name" content="{TITLE}">
<meta property="og:title" content="{TITLE}">
<meta property="og:description" content="{DESC}">
<meta property="og:url" content="{SITE}/">
<meta property="og:image" content="{IMAGE}">
<meta name="twitter:card" content="summary">
<meta name="twitter:title" content="{TITLE}">
<meta name="twitter:description" content="{DESC}">
<meta name="twitter:image" content="{IMAGE}">
<meta name="description" content="{DESC}">
"""


def main():
    index = Path(streamlit.__file__).parent / "static" / "index.html"
    html = index.read_text()
    if MARK in html:
        print(f"already patched: {index}")
        return
    html = html.replace("<head>", "<head>" + TAGS, 1)
    # crawlers see the pre-JS title; page_config only fixes it client-side
    html = html.replace("<title>Streamlit</title>", f"<title>{TITLE}</title>", 1)
    index.write_text(html)
    print(f"patched: {index}")


if __name__ == "__main__":
    main()

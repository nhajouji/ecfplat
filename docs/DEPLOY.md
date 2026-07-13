# Deploying to elliptic-curves.info

Streamlit Community Cloud cannot serve a custom domain, so the site moves to
**Render** (the repo ships a `render.yaml` blueprint). Query-param routing
works identically there — nothing in the app changes.

## One-time setup (~15 minutes, needs your accounts)

### 1. Create the Render service
1. Sign in at https://render.com (GitHub login is easiest).
2. **New + → Blueprint**, pick the `nhajouji/ecfplat` repo.
   Render reads `render.yaml`: Python 3.13.5, `pip install -r
   requirements.txt`, and the streamlit start command on `$PORT`.
3. Plan: the blueprint says `starter` (~$7/mo, always on). The free plan also
   works but the app spins down after 15 idle minutes and cold-starts in
   ~30-60 s — fine for testing, annoying for visitors.
4. Deploy; check the `*.onrender.com` URL works (drill into
   `/Explorer?a=6&p=101` to confirm query-param routing survives).

### 2. Point elliptic-curves.info at it
1. In the Render service: **Settings → Custom Domains → Add**, enter both
   `elliptic-curves.info` and `www.elliptic-curves.info`.
2. Render shows the DNS records it wants. At your registrar's DNS panel add:
   - **Apex** (`elliptic-curves.info`): `A` record → `216.24.57.1`
     (Render's apex IP — use whatever value the dashboard shows), or an
     `ALIAS`/`ANAME` record to the `.onrender.com` host if the registrar
     supports it.
   - **www**: `CNAME` → `<service>.onrender.com`.
3. Wait for DNS to propagate (minutes to an hour). Render provisions the TLS
   certificate automatically once it verifies the records; both hosts then
   serve over https.

### 3. Retire the old deployment
Delete (or just stop sharing) the Streamlit Community Cloud app so there is a
single canonical URL.

## Notes
- **Python version** is pinned twice: `.python-version` (3.13) and the
  `PYTHON_VERSION` env var in `render.yaml`. Keep them in sync; 3.14 is
  untested for this stack.
- **Memory**: the app eagerly loads only the small stores (~2 MB); the 83 MB
  `ecqf_ord_pcbij_ext.json` is *not* loaded by the site. 512 MB instances are
  plenty.
- **requirements.txt is deploy-only**: matplotlib/plotly were dropped because
  no page imports them anymore. Notebooks and dev tools that want them can
  install them locally.
- **Updating the site** stays git-push-driven: Render auto-deploys on every
  push to the connected branch.

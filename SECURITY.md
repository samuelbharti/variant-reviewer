# Security Policy

## Supported versions

I maintain this project on my own. I fix security problems on `main` and in the
next release. I do not backport fixes to older tags.

## Reporting a problem

**Please do not open a public issue for a security problem.**

Email me at <samuelbharti.io@gmail.com>. Tell me what you found and, if you
can, how to reproduce it. I will acknowledge your report within a few days and
tell you what I plan to do about it.

## Worth knowing before you report

- The app queries public, read-only bioinformatics APIs. None of them needs a
  key, and the app sends no user data to them beyond the gene or variant you
  search for.
- The AI assistant is bring your own key. A key you paste stays in server
  memory for your session only. The app never writes it to disk or to a log.
- If you deploy the app yourself, the keys and the environment you deploy into
  are yours to secure. Serve it over HTTPS, since a pasted key travels from the
  browser to the server.

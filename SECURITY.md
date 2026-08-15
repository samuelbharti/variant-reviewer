# Security Policy

## Supported versions

Fixes land on `main` and go out with the next release. Older tags are not
patched.

## Reporting a problem

Please do not open a public issue for a security problem. Email
<samuelbharti.io@gmail.com> instead, describing what you found and, where you
can, the steps to reproduce it. You will get an acknowledgement within a few
days, along with what happens next.

## Worth knowing before you report

- The app queries public, read-only bioinformatics APIs. None of them needs a
  key, and the app sends no user data to them beyond the gene or variant you
  search for.
- The AI assistant is bring your own key. A key you paste stays in server
  memory for your session only. The app never writes it to disk or to a log.
- If you deploy the app yourself, the keys and the environment you deploy into
  are yours to secure. Serve it over HTTPS, since a pasted key travels from the
  browser to the server.

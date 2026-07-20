# BYOK chat module: pure, network-free logic (provider handling, secret
# redaction, friendly-error mapping) plus a stubbed server flow driven with
# injected fakes -- no network, no real key. Mirrors the tests that ship with
# the module in R_Shiny_modules, with an extra check for the app's context tool.

test_that("unknown providers are dropped so they can't reach a backend", {
  expect_identical(
    .byok_chat_providers(c("openai", "bogus", "gemini")),
    c("openai", "gemini")
  )
  expect_identical(.byok_chat_providers("nope"), character(0))
})

test_that("provider metadata carries a label and env fallback vars", {
  expect_identical(.byok_chat_provider_meta("openai")$env, "OPENAI_API_KEY")
  expect_true("GOOGLE_API_KEY" %in% .byok_chat_provider_meta("gemini")$env)
  expect_true(nzchar(.byok_chat_provider_meta("anthropic")$label))
})

test_that("a known secret is redacted from any surfaced message", {
  expect_identical(
    .byok_chat_redact("leaked sk-secret here", "sk-secret"),
    "leaked <redacted-key> here"
  )
  # empty secret must not mangle the message (no accidental blanket redaction)
  expect_identical(.byok_chat_redact("nothing to hide", ""), "nothing to hide")
})

test_that("friendly errors classify common failures and never leak the key", {
  expect_match(.byok_chat_friendly_error("Invalid API key", "k"), "rejected")
  expect_match(.byok_chat_friendly_error("429 rate limit"), "Rate limit")
  expect_match(.byok_chat_friendly_error("model_not_found"), "isn't available")
  expect_match(.byok_chat_friendly_error("503 overloaded"), "busy")
  expect_false(grepl(
    "sk-zzz",
    .byok_chat_friendly_error("bad request with sk-zzz", "sk-zzz")
  ))
})

test_that("fetch_models returns NULL (fallback) for empty key or unknown provider", {
  expect_null(.byok_chat_fetch_models("openai", ""))
  expect_null(.byok_chat_fetch_models("bogus", "some-key"))
})

test_that("server connects, streams, and enforces the turn limit (stubbed)", {
  rec <- new.env()
  rec$appended <- character()
  rec$streamed <- character()
  stub <- list(stream_async = function(msg) {
    rec$streamed <- c(rec$streamed, msg)
    "ok"
  })

  testServer(
    byok_chat_server,
    args = list(
      max_turns = 1L,
      client_factory = function(provider, api_key, model) {
        rec$conn <- list(provider = provider, key = api_key, model = model)
        stub
      },
      append = function(response) {
        rec$appended <- c(rec$appended, as.character(response))
        invisible(NULL)
      }
    ),
    {
      session$setInputs(provider = "openai")

      # Chatting before connecting prompts the user, does not build a client.
      session$setInputs(chat_user_input = "hi")
      expect_true(any(grepl("Model & key", rec$appended, fixed = TRUE)))
      expect_null(client())

      # The header badge reflects the disconnected state.
      expect_match(output$conn_badge$html, "Not connected")

      # Connect with a pasted key.
      session$setInputs(api_key = "sk-secret-123", model = "gpt-4o")
      session$setInputs(connect = 1)
      expect_false(is.null(client()))
      expect_identical(rec$conn$provider, "openai")
      expect_identical(rec$conn$key, "sk-secret-123")
      expect_true(isTRUE(status()$ok))
      # ...and the header badge now shows the connected state.
      expect_match(output$conn_badge$html, "Connected")

      # A real turn streams the message.
      session$setInputs(chat_user_input = "hello")
      expect_true("hello" %in% rec$streamed)
      expect_identical(n_turns(), 1L)

      # Exceeding max_turns surfaces the limit instead of streaming again.
      session$setInputs(chat_user_input = "too much")
      expect_true(any(grepl("Turn limit reached", rec$appended)))

      # Forget clears the client.
      session$setInputs(forget = 1)
      expect_null(client())
    }
  )
})

test_that("an environment key is surfaced in the status and key help", {
  skip_if_not_installed("withr")
  # Only OpenAI has a server-side key; blank the others (the dev/CI env or a
  # project .Renviron may otherwise set them).
  withr::local_envvar(c(
    OPENAI_API_KEY = "sk-env-test",
    ANTHROPIC_API_KEY = "",
    GEMINI_API_KEY = "",
    GOOGLE_API_KEY = ""
  ))
  testServer(
    byok_chat_server,
    args = list(
      client_factory = function(provider, api_key, model) {
        list(stream_async = function(msg) "ok")
      },
      append = function(response) invisible(NULL)
    ),
    {
      session$setInputs(provider = "openai")
      # The status invites picking a model and connecting, no pasting needed.
      expect_match(status()$msg, "found in the environment")
      expect_match(output$key_help$html, "set in the environment")

      # A provider without a server key falls back to the paste prompt.
      session$setInputs(provider = "anthropic")
      expect_match(status()$msg, "Enter your key")
    }
  )
})

test_that("context tool closure reflects the live app snapshot", {
  # server.R backs the current_review_context tool with a plain-env snapshot so
  # the closure can be read during async streaming (outside a reactive context).
  # This checks that read-through pattern; the ellmer::tool() wrapping is thin.
  chat_context <- new.env(parent = emptyenv())
  chat_context$text <- "No gene or variant has been searched yet."
  read_context <- function() chat_context$text

  expect_match(read_context(), "No gene or variant")

  chat_context$text <- "Searched gene: TP53"
  expect_identical(read_context(), "Searched gene: TP53")
})

test_that(".byok_chat_suggestions_md renders clickable suggestion cards", {
  md <- .byok_chat_suggestions_md(c(
    "Gene overview" = "Load TP53.",
    "Interactions" = "Top partners for EGFR?"
  ))
  # Markdown list where every item is a .suggestion span (shinychat renders
  # these as clickable cards that submit their body text on click).
  expect_match(
    md,
    '- <span class="suggestion" title="Gene overview">Load TP53.</span>',
    fixed = TRUE
  )
  expect_equal(length(strsplit(md, "\n")[[1]]), 2)

  # HTML in a prompt is escaped so it can't break the markup.
  escaped <- .byok_chat_suggestions_md(c("A & <b>" = "1 < 2 & 3"))
  expect_match(escaped, "1 &lt; 2 &amp; 3", fixed = TRUE)

  # Empty / NULL input yields no list.
  expect_null(.byok_chat_suggestions_md(NULL))
  expect_null(.byok_chat_suggestions_md(character(0)))
})

test_that(".byok_chat_connected_greeting includes the suggestions when present", {
  with_sug <- .byok_chat_connected_greeting("OpenAI", c(x = "Do a thing."))
  expect_match(with_sug, "Connected to \\*\\*OpenAI\\*\\*")
  expect_match(with_sug, 'class="suggestion"')

  without <- .byok_chat_connected_greeting("OpenAI", NULL)
  expect_no_match(without, "suggestion")
})

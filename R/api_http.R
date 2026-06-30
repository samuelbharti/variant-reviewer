# Shared HTTP helper for all external API clients.
#
# Every client funnels through vr_api_get() so timeouts, retries, and error
# handling live in one place. The return value is always a normalized list so
# callers never need their own tryCatch():
#   list(ok = TRUE,  status = 200L, data = <parsed JSON>)
#   list(ok = FALSE, status = <int|NA>, error = "<message>", data = NULL)

# A friendly label for the source of a failed request, used in error messages.
vr_api_get <- function(
  base_url,
  path = NULL,
  query = list(),
  source = "API",
  timeout = 15,
  max_tries = 3
) {
  req <- httr2::request(base_url)
  if (!is.null(path)) {
    req <- httr2::req_url_path_append(req, path)
  }
  if (length(query) > 0) {
    # Drop NULL/empty query values so we never send "param=".
    query <- query[!vapply(query, is_blank, logical(1))]
    req <- do.call(httr2::req_url_query, c(list(req), query))
  }
  req <- httr2::req_timeout(req, timeout)
  req <- httr2::req_retry(req, max_tries = max_tries)
  # Don't let httr2 raise on HTTP errors; we normalize them ourselves.
  req <- httr2::req_error(req, is_error = function(resp) FALSE)
  req <- httr2::req_user_agent(req, "variant-reviewer (Shiny app)")

  tryCatch(
    {
      resp <- httr2::req_perform(req)
      status <- httr2::resp_status(resp)
      if (status >= 200 && status < 300) {
        list(
          ok = TRUE,
          status = status,
          # check_type = FALSE: some sources (e.g. STRING) send "text/json".
          data = httr2::resp_body_json(
            resp,
            check_type = FALSE,
            simplifyVector = FALSE
          ),
          error = NULL
        )
      } else {
        list(
          ok = FALSE,
          status = status,
          data = NULL,
          error = paste0(source, " returned HTTP ", status, ".")
        )
      }
    },
    error = function(e) {
      list(
        ok = FALSE,
        status = NA_integer_,
        data = NULL,
        error = paste0("Could not reach ", source, ": ", conditionMessage(e))
      )
    }
  )
}

# TRUE for NULL, NA, empty string, or zero-length values.
is_blank <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    (length(x) == 1 && (is.na(x) || identical(trimws(as.character(x)), "")))
}

# Pull a value from a nested list by key path, returning `default` if any level
# is missing or NULL. e.g. pluck_at(x, "ensembl", "gene").
pluck_at <- function(x, ..., default = NULL) {
  keys <- c(...)
  for (key in keys) {
    if (is.null(x) || !is.list(x) || is.null(x[[key]])) {
      return(default)
    }
    x <- x[[key]]
  }
  if (is.null(x)) default else x
}

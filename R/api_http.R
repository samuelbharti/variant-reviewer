# Shared HTTP helpers for all external API clients.
#
# Every client funnels through vr_api_get() (REST) or vr_api_post_json()
# (GraphQL/JSON POST) so timeouts, retries, caching, and error handling live in
# one place. The return value is always a normalized list so callers never need
# their own tryCatch():
#   list(ok = TRUE,  status = 200L, data = <parsed JSON>)
#   list(ok = FALSE, status = <int|NA>, error = "<message>", data = NULL)
#
# Successful responses are cached in-process for VR_CACHE_TTL seconds so
# repeated searches for the same gene/variant are instant; failures are never
# cached, so a transient outage doesn't get stuck.

VR_CACHE_TTL <- 1800 # 30 minutes
vr_cache <- cachem::cache_mem(max_age = VR_CACHE_TTL)

# Stable cache key for a request (lowercase hex hash -> valid cachem key).
vr_cache_key <- function(...) {
  rlang::hash(list(...))
}

# GET a REST endpoint. `source` is a friendly label used in error messages.
vr_api_get <- function(
  base_url,
  path = NULL,
  query = list(),
  source = "API",
  timeout = 15,
  max_tries = 3
) {
  if (length(query) > 0) {
    # Drop NULL/empty query values so we never send "param=".
    query <- query[!vapply(query, is_blank, logical(1))]
  }
  key <- vr_cache_key("GET", base_url, path, query)

  vr_cached(key, function() {
    req <- httr2::request(base_url)
    if (!is.null(path)) {
      req <- httr2::req_url_path_append(req, path)
    }
    if (length(query) > 0) {
      req <- do.call(httr2::req_url_query, c(list(req), query))
    }
    req <- vr_req_defaults(req, timeout, max_tries)
    vr_perform(req, source)
  })
}

# POST a JSON body (e.g. a GraphQL query) and parse the JSON response.
vr_api_post_json <- function(
  url,
  body,
  source = "API",
  timeout = 20,
  max_tries = 3
) {
  key <- vr_cache_key("POST", url, body)

  vr_cached(key, function() {
    req <- httr2::request(url)
    req <- httr2::req_body_json(req, body)
    req <- vr_req_defaults(req, timeout, max_tries)
    vr_perform(req, source)
  })
}

# Return a cached successful result for `key`, otherwise run `fetch()` and cache
# it only when it succeeded.
vr_cached <- function(key, fetch) {
  hit <- vr_cache$get(key, missing = NULL)
  if (!is.null(hit)) {
    return(hit)
  }
  res <- fetch()
  if (isTRUE(res$ok)) {
    vr_cache$set(key, res)
  }
  res
}

# Apply the shared request options (timeout, retries, no-raise, user agent).
vr_req_defaults <- function(req, timeout, max_tries) {
  req <- httr2::req_timeout(req, timeout)
  req <- httr2::req_retry(req, max_tries = max_tries)
  # Don't let httr2 raise on HTTP errors; we normalize them ourselves.
  req <- httr2::req_error(req, is_error = function(resp) FALSE)
  httr2::req_user_agent(req, "variant-reviewer (Shiny app)")
}

# Perform a request and normalize the result into the standard ok/status/data
# /error list.
vr_perform <- function(req, source) {
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

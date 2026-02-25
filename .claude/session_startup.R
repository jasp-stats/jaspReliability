# JASP Module - R Session Startup for Claude Code
# Run this script in your interactive R session (RStudio/Positron/radian)
# to prepare and hand over the session to Claude Code.
#
# Usage: source(".claude/session_startup.R")

# Fix cli::get_spinner() conflict with testthat in btw/evaluate context
options(cli.spinner = "line")

# Fix locale issue with renv.lock files created on non-English systems
if (.Platform$OS.type == "windows") {
  Sys.setlocale("LC_ALL", "English_United States.utf8")
} else {
  Sys.setlocale("LC_ALL", "C.UTF-8")
}

# Install order matches container_entrypoint.sh:
# 1. Install injected packages (btw, mcptools) first
# 2. renv::restore() so lockfile-pinned versions win for shared deps
# 3. Install the module + jaspTools (already handled by restore for deps)
renv::install(c("btw", "mcptools"), prompt = FALSE)
renv::restore(prompt = FALSE)
renv::install(c(".", "jasp-stats/jaspTools"), prompt = FALSE)
library(jaspTools)
setupJaspTools()
setPkgOption("module.dirs", ".")
setPkgOption("reinstall.modules", FALSE)
btw::btw_mcp_session()

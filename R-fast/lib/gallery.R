# A single HTML page showing every tax figure: one section per area (United States first, then
# states alphabetically) with its four figures, and a left sidebar linking to each. Each area's
# figures sit in 1, 2 or 4 columns by window width (phone, laptop, wide monitor), so an area is
# never split 3 + 1; on narrow windows the sidebar becomes a bar of abbreviations at the top.
# Images are the SVG files, referenced by relative path and loaded lazily, so the page opens fast
# and stays sharp at any zoom.
#
# figs: data frame from save_tax_plots() with columns fips, abbr, figure and svg.
# fig_rel_dir: the SVG folder relative to the page. subtitle: a line describing the run.

write_tax_gallery <- function(path, figs, fig_rel_dir, subtitle) {
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  labels <- c(LYG = "Life-years gained (cumulative)", DeathAvoided = "Deaths avoided (cumulative)",
              Prev = "Smoking prevalence", PrevRed = "Prevalence reduction")
  full <- c(setNames(usmap::statepop$full, usmap::statepop$abbr), US = "United States")
  figs$name <- unname(full[figs$abbr])
  areas <- unique(figs[order(figs$abbr != "US", figs$name), c("abbr", "name")])

  nav <- paste0('<a href="#', areas$abbr, '"><span class="ab">', areas$abbr, '</span><span class="nm">',
                esc(areas$name), "</span></a>", collapse = "\n      ")
  sections <- vapply(seq_len(nrow(areas)), function(i) {
    st <- figs[figs$abbr == areas$abbr[i], ]
    st <- st[order(match(st$figure, names(labels))), ]
    cards <- paste0('<figure><img src="', esc(file.path(fig_rel_dir, st$svg)), '" alt="', esc(areas$name[i]), " - ",
                    esc(labels[st$figure]), '" loading="lazy" width="600" height="400"><figcaption>',
                    esc(labels[st$figure]), "</figcaption></figure>", collapse = "\n        ")
    sprintf('    <section id="%s">\n      <h2>%s <span>%s</span></h2>\n      <div class="grid">\n        %s\n      </div>\n    </section>',
            areas$abbr[i], esc(areas$name[i]), areas$abbr[i], cards)
  }, character(1))

  html <- c('<!doctype html>', '<html lang="en">', '<head>', '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            '<title>Tax policy figures</title>', '<style>',
            ':root { --bg: #f6f6f4; --card: #ffffff; --text: #1d1d1f; --muted: #6b6b68; --border: #e2e2de; --hover: #e8e8e4; --side: 200px; }',
            '* { box-sizing: border-box; }',
            'body { margin: 0; background: var(--bg); color: var(--text); font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }',
            'nav { position: fixed; top: 0; left: 0; bottom: 0; width: var(--side); overflow-y: auto; padding: 24px 12px 24px 16px; border-right: 1px solid var(--border); background: var(--bg); }',
            'nav a { display: flex; gap: 8px; padding: 3px 8px; border-radius: 4px; color: var(--text); text-decoration: none; font-size: 13px; }',
            'nav a:hover { background: var(--hover); }',
            'nav .ab { width: 22px; color: var(--muted); font-size: 12px; }',
            'nav a:first-child { margin-bottom: 8px; font-weight: 600; }',
            '.content { margin-left: var(--side); padding: 28px 32px 64px; max-width: calc(2600px + var(--side)); }',
            'h1 { font-size: 22px; font-weight: 600; margin: 0 0 4px; }',
            'header p { margin: 0; color: var(--muted); font-size: 14px; }',
            'section { margin-top: 32px; scroll-margin-top: 16px; }',
            'h2 { font-size: 18px; font-weight: 600; margin: 0 0 12px; }',
            'h2 span { color: var(--muted); font-weight: 400; font-size: 14px; margin-left: 6px; }',
            '.grid { display: grid; grid-template-columns: 1fr; gap: 16px; }',
            '@media (min-width: 1200px) { .grid { grid-template-columns: repeat(2, 1fr); } }',
            '@media (min-width: 2200px) { .grid { grid-template-columns: repeat(4, 1fr); } }',
            'figure { margin: 0; background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 12px 12px 8px; }',
            'figure img { display: block; width: 100%; height: auto; }',
            'figcaption { font-size: 13px; color: var(--muted); margin-top: 4px; }',
            '@media (max-width: 900px) {',
            '  nav { position: sticky; bottom: auto; width: auto; z-index: 1; display: flex; flex-wrap: wrap; gap: 2px; padding: 8px 16px; border-right: 0; border-bottom: 1px solid var(--border); background: rgba(246, 246, 244, 0.95); }',
            '  nav a { padding: 2px 6px; } nav a:first-child { margin-bottom: 0; } nav .nm { display: none; } nav .ab { width: auto; color: var(--text); }',
            '  .content { margin-left: 0; padding: 20px 16px 48px; } section { scroll-margin-top: 110px; }',
            '}',
            '</style>', '</head>', '<body>',
            '  <nav>', paste0("      ", nav), '  </nav>',
            '  <div class="content">',
            '    <header>', '      <h1>Cigarette tax policy: figures</h1>', paste0("      <p>", esc(subtitle), "</p>"), '    </header>',
            sections, '  </div>', '</body>', '</html>')
  writeLines(html, path, useBytes = TRUE)
  invisible(path)
}

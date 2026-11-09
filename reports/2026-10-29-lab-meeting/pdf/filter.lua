-- Pandoc filter for REPORT.pdf (see build.py).
--  * images: paths made relative to pdf/, sized from their pixel width, height capped
--  * horizontal rules -> page breaks; \FloatBarrier before every heading
--  * blockquotes -> boxed notes in the template's own mdframed styles

local sizes = {}
do
  local f = io.open('build/imgsizes.json', 'r')
  if f then sizes = pandoc.json.decode(f:read('a')); f:close() end
end
local LINEWIDTH_IN = 6.9   -- letter paper, 0.8 in margins (annotation.sty)
local MIN_DPI = 90         -- small gel scans are not blown up past this

function Image(img)
  local src = img.src
  if src:match('%.svg$') then
    img.src = 'spring_types.pdf'
    img.attributes.width = '70%'
    return img
  end
  local s = sizes[src]
  img.src = '../' .. src
  if s then
    local w = math.min(0.85 * LINEWIDTH_IN, s[1] / MIN_DPI)
    img.attributes.width = string.format('%.2fin', w)
    img.attributes.height = '38%'
  end
  return img
end

-- tables: wide ones get wrapping columns with widths proportional to their longest cell
local function cell_len(cell)
  local s = pandoc.utils.stringify(cell.contents)
  return utf8.len(s) or #s
end

local function longest_word(cell)
  local m = 0
  for w in pandoc.utils.stringify(cell.contents):gmatch('%S+') do m = math.max(m, utf8.len(w) or #w) end
  return m
end

local CHARS_PER_LINE = 100   -- \small cmbright across the text width, roughly

function Table(tbl)
  local ncol = #tbl.colspecs
  local len, raw, word = {}, {}, {}
  for c = 1, ncol do len[c] = 3; raw[c] = 0; word[c] = 0 end
  local function scan(rows)
    for _, row in ipairs(rows) do
      for c, cell in ipairs(row.cells) do
        if c <= ncol then
          local n = cell_len(cell)
          raw[c] = math.max(raw[c], n)
          len[c] = math.max(len[c], math.min(n, 90))
          word[c] = math.max(word[c], longest_word(cell))
        end
      end
    end
  end
  scan(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do scan(body.body) end
  local total, rawtotal = 0, 0
  for c = 1, ncol do total = total + len[c]; rawtotal = rawtotal + raw[c] end
  if rawtotal <= 75 then return nil end    -- every row fits on one line: natural column widths
  -- proportional to content, but never narrower than the column's longest word
  local f, sum = {}, 0
  for c = 1, ncol do
    f[c] = math.max(len[c] / total, (word[c] + 3) / CHARS_PER_LINE)
    sum = sum + f[c]
  end
  for c = 1, ncol do
    tbl.colspecs[c] = { tbl.colspecs[c][1], 0.98 * f[c] / sum }
  end
  return tbl
end

function HorizontalRule()
  return pandoc.RawBlock('latex', '\\newpage')
end

function Header(h)
  return { pandoc.RawBlock('latex', '\\FloatBarrier'), h }
end

local function drop_lead(inl)
  table.remove(inl, 1)
  if inl[1] and inl[1].t == 'Space' then table.remove(inl, 1) end
end

function BlockQuote(bq)
  local title, style = 'Note', 'mdorangebox'
  local first = bq.content[1]
  if first and first.t == 'Para' then
    local inl = first.content
    if inl[1] and inl[1].t == 'Strong' then
      title = pandoc.utils.stringify(inl[1]):gsub('%.$', '')
      drop_lead(inl)
    elseif inl[1] and inl[1].t == 'Str' and inl[1].text == 'Note:' then
      drop_lead(inl)
    end
  end
  if title:match('^Milestone') then style = 'mdgreenbox'
  elseif title:match('^How to reuse') then style = 'mdbluebox'
  elseif title:match('^Methods detail') then style = 'mdgreybox' end
  local out = { pandoc.RawBlock('latex', '\\begin{reportbox}{' .. style .. '}{' .. title .. '}') }
  for _, b in ipairs(bq.content) do table.insert(out, b) end
  table.insert(out, pandoc.RawBlock('latex', '\\end{reportbox}'))
  return out
end

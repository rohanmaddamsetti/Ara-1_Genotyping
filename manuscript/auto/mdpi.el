(TeX-add-style-hook "mdpi"
 (lambda ()
    (LaTeX-add-environments
     "thebibliography")
    (TeX-add-symbols
     '("corres" ["argument"] 1)
     '("address" ["argument"] 1)
     '("heading" 1)
     '("keyword" 1)
     '("MSC" 1)
     '("PACS" 1)
     '("history" 1)
     '("pubyear" 1)
     '("pubvolume" 1)
     '("doinum" 1)
     '("lastpage" 1)
     '("Author" 1)
     '("Title" 1)
     '("cor" 2)
     "org"
     "firstargument"
     "corresfirstargument"
     "journalname"
     "logosize"
     "logo"
     "issn"
     "web"
     "logolength"
     "header"
     "abstractkeywords"
     "cright"
     "citeyearpar"
     "natexlab"
     "maketitle"
     "sectionmark"
     "subsectionmark"
     "makelabel")
    (TeX-run-style-hooks
     "url"
     "hyperref"
     "natbib"
     "sort&compress"
     "color"
     "caption"
     "float"
     "lineno"
     "ifthen"
     "lastpage"
     "soul"
     "graphicx"
     "fancyhdr"
     "times"
     "indentfirst"
     "calc"
     "art12"
     "article"
     "12pt")))


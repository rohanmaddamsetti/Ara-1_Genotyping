(TeX-add-style-hook "genotyping"
 (lambda ()
    (LaTeX-add-bibitems
     "ref-journal"
     "ref-book")
    (TeX-run-style-hooks
     "latex2e"
     "mdpi12"
     "mdpi"
     "genes"
     "article"
     "submit"
     "moreauthors"
     "pdftex"
     "12pt"
     "a4paper")))


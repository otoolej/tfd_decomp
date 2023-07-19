(TeX-add-style-hook
 "ta"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "booktabs")
   (LaTeX-add-labels
    "tab:combined"))
 :latex)


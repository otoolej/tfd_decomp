
(defun org-table-transpose ()
  "Transpose the Org table at point."
  (interactive)
  (when (org-at-table-p)
    (let* ((table (org-table-to-lisp))
           (transposed (apply #'cl-mapcar #'list table)))
      (delete-region (org-table-begin) (org-table-end))
      (insert (mapconcat (lambda (row)
                           (concat "| " (mapconcat 'identity row " | ") " |"))
                         transposed
                         "\n"))
      (org-table-align))))

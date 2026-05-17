csee_talk.html: csee_talk.qmd
	quarto render $<

## retrieve revealjs-edited version
switch:
	mv -b ~/Downloads/csee_talk.qmd .

quarto_install:
	quarto add emilhvitfeldt/quarto-revealjs-editable

upload: csee_talk.html
	scp csee_talk.html ms.mcmaster.ca:~/public_html/misc/pv_csee.html

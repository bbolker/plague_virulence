ssc_talk.html: ssc_talk.qmd
	quarto render $<

csee_talk.html: csee_talk.qmd
	quarto render $<

## retrieve revealjs-edited version
switch:
	mv -b ~/Downloads/ssc_talk.qmd .

quarto_install:
	quarto add emilhvitfeldt/quarto-revealjs-editable

upload: ssc_talk.html
	scp ssc_talk.html ms.mcmaster.ca:~/public_html/misc/ssc_plaguetalk_2026.html

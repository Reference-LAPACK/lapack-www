

all: bug_list.html contributor-list.html faq.html release_notes.html improvement.html index.html err coding lawn

bug_list.html: bug_list.txt
	asciidoc -a toc bug_list.txt

contributor-list.html: contributor-list.txt
	asciidoc  -a toc contributor-list.txt

faq.html: faq.txt
	asciidoc -a toc faq.txt

release_notes.html: release_notes.txt
	asciidoc release_notes.txt

improvement.html: improvement.txt
	asciidoc improvement.txt

index.html: index.txt
	asciidoc -a toc -a toc-title="Menu" index.txt

err: Errata/index2.txt
	@(cd Errata && make && cd ..)

lawn: lawns/index.txt
	@(cd lawns && make && cd ..)

coding: lapack-coding/program-style.txt
	@(cd lapack-coding && make && cd ..)

publish:
	scp *.txt *.html netlib.org:/netlib/lapack
	scp Errata/*.txt Errata/*.html netlib.org:/netlib/lapack/Errata
	scp lawns/*.txt lawns/*.html lawns/lawn.bib netlib.org:/netlib/lapack/lawns
	scp lapack-coding/*.txt lapack-coding/*.html netlib.org:/netlib/lapack-dev/lapack-coding

pub_bug:
	scp bug_list.* netlib.org:/netlib/lapack
	scp Errata/*.txt Errata/*.html netlib.org:/netlib/lapack/Errata

clean:
	rm -rf *.html Errata/*.html lapack-coding/*.html
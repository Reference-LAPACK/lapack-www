

all: bug_list.html contributor-list.html faq.html release_notes improvement.html index.html err coding

bug_list.html: bug_list.txt
	asciidoc -a toc bug_list.txt

contributor-list.html: contributor-list.txt
	asciidoc  -a toc contributor-list.txt

faq.html: faq.txt
	asciidoc -a toc faq.txt

release_notes: release_notes.txt release_notes-3.3.0.txt lapack-3.3.1.txt lapack-3.3.0.txt lapack-3.4.0.txt lapacke.txt errata_from_331_to_340.txt
	asciidoc release_notes.txt
	asciidoc release_notes-3.3.0.txt
	asciidoc lapack-3.3.1.txt
	asciidoc lapack-3.3.0.txt
	asciidoc -a toc lapack-3.4.0.txt
	asciidoc -a toc lapacke.txt
	asciidoc -a toc errata_from_331_to_340.txt

improvement.html: improvement.txt
	asciidoc improvement.txt

index.html: index.txt
	asciidoc -a toc -a toc-title="Menu" index.txt

err: Errata/index.txt Errata/errata_scalapack.txt
	@(cd Errata && make && cd ..)

lawn: lawns/index.txt
	@(cd lawns && make && cd ..)
	scp lawns/*.txt lawns/*.html lawns/lawn.bib netlib.org:/netlib/lapack/lawns

coding: lapack-coding/program-style.txt
	@(cd lapack-coding && make && cd ..)

publish:
	scp *.txt *.html netlib.org:/netlib/lapack
	scp Errata/*.txt Errata/*.html netlib.org:/netlib/lapack/Errata
	scp lapack-coding/*.txt lapack-coding/*.html netlib.org:/netlib/lapack-dev/lapack-coding

pub_bug: bug_list.html err
	scp bug_list.* netlib.org:/netlib/lapack
	scp Errata/*.txt Errata/*.html netlib.org:/netlib/lapack/Errata
	scp -r Errata/vrac/* netlib.org:/netlib/lapack/Errata/vrac

clean:
	rm -rf *.html Errata/*.html lapack-coding/*.html

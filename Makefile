

all: bug_list.html contributor-list.html faq.html release_notes improvement.html index err coding

bug_list.html: bug_list.txt
	asciidoc -a toc bug_list.txt

contributor-list.html: contributor-list.txt
	asciidoc  -a toc contributor-list.txt

faq.html: faq.txt
	asciidoc -a toc faq.txt

release_notes: release_notes.txt release_notes-3.9.0.txt lapack-3.9.0.txt

	asciidoc release_notes-3.9.0.txt
	asciidoc -a toc lapack-3.9.0.txt
	
	scp release_notes-3.9.0.html release_notes-3.9.0.txt julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack
	scp lapack-3.9.0.html lapack-3.9.0.txt julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack


improvement.html: improvement.txt
	asciidoc improvement.txt

index: index.txt
	asciidoc -a toc -a toc-title="Menu" index.txt
	scp index.html index.txt julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack
	

err: Errata/index2.txt Errata/index2.txt Errata/errata_scalapack.txt
	@(cd Errata && make && cd ..)

lawn: lawns/index.txt
	@(cd lawns && make && cd ..)
	scp lawns/*.txt lawns/*.html lawns/lawn.bib julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack/lawns

coding: lapack-coding/program-style.txt
	@(cd lapack-coding && make && cd ..)

publish:
	scp *.txt *.html zoot.icl.utk.edu:/nfs/www/netlib/lapack
	scp Errata/*.txt Errata/*.html julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack/Errata
	scp lapack-coding/*.txt lapack-coding/*.html julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack-dev/lapack-coding

pub_bug: bug_list.html err
	scp bug_list.* zoot.icl.utk.edu:/nfs/www/netlib/lapack
	scp Errata/*.txt Errata/*.html julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack/Errata
	scp -r Errata/vrac/*  julie@zoot.icl.utk.edu:/nfs/www/netlib/lapack/Errata/vrac

clean:
	rm -rf *.html Errata/*.html lapack-coding/*.html

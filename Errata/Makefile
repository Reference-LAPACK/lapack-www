
all: index2.html errata_scalapack.html

index2.html: index2.txt
	asciidoc -a toc -a toc-title="Menu" index2.txt
	asciidoc -a toc -a toc-title="Menu" errata_3.3.0.txt

errata_scalapack.html: errata_scalapack.txt
	asciidoc -a toc -a toc-title="Menu" errata_scalapack.txt

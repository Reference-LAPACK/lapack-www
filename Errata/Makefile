
all: index.html errata_scalapack.html

index.html: index.txt
	asciidoc -a toc -a toc-title="Menu" index.txt
	asciidoc -a toc -a toc-title="Menu" errata_3.3.0.txt

errata_scalapack.html: errata_scalapack.txt
	asciidoc -a toc -a toc-title="Menu" errata_scalapack.txt

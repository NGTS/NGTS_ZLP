all: process.html schematic.png schematic.svg schematic.pdf

%.png: %.dot
	dot -T png $< -o $@

%.svg: %.dot
	dot -T svg $< -o $@

%.pdf: %.dot
	dot -T pdf $< -o $@

%.html: %.asciidoc
	asciidoctor $<

.PHONY: clean

clean:
	@rm -f *.pdf *.png *.svg

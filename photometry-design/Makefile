all: schematic.png schematic.svg schematic.pdf

%.png: %.dot
	dot -T png $< -o $@

%.svg: %.dot
	dot -T svg $< -o $@

%.pdf: %.dot
	dot -T pdf $< -o $@

#----------------------------------------------------------------------
# Build the Tacode version 2 report.
#
# The English version is the authoritative one; the Japanese version follows it.
#
#   make            both PDFs
#   make en         the English one only
#   make ja         the Japanese one only
#   make figure     rebuild every figure from the data in the repository
#   make check      check the sources without typesetting (no TeX needed)
#   make clean      remove the intermediate files
#
# The English version builds with pdflatex, the Japanese one with platex and
# dvipdfmx, which is the toolchain the version 1 manuscript used. Both need the
# figures under figure/; they are committed, so a plain `make` needs no Python.
#
# The Japanese fonts have to be selected once before they are embedded, or the
# PDF carries only font names and shows nothing on a machine without them:
#
#   kanji-config-updmap-user haranoaji        # this user
#   sudo kanji-config-updmap-sys haranoaji    # everyone on the machine
#
# `pdffonts tacode_v2_ja.pdf` should then show "emb yes" for every font.
#----------------------------------------------------------------------

PYTHON  ?= python3
PDFLATEX?= pdflatex
PLATEX  ?= platex
DVIPDF  ?= dvipdfmx -V 7

.PHONY: all en ja figure check clean

all: en ja

en: tacode_v2_en.pdf
ja: tacode_v2_ja.pdf

tacode_v2_en.pdf: tacode_v2_en.tex
	$(PDFLATEX) $<
	$(PDFLATEX) $<

tacode_v2_ja.pdf: tacode_v2_ja.tex
	$(PLATEX) $<
	$(PLATEX) $<
	$(DVIPDF) $(basename $<)

figure:
	$(PYTHON) generate_figure.py

check:
	$(PYTHON) check_tex.py

clean:
	rm -f *.aux *.log *.out *.dvi *.toc

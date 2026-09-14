#!/usr/bin/env python3
#
# レポートの .tex を**組版せずに**検査する。
#
# TeX 処理系が入っていない環境でも、typo の大半はここで捕まえられる:
# 環境の対応、括弧と $ の数、\ref に対応する \label、\cite に対応する \bibitem、
# \includegraphics が指すファイルの実在、そして Markdown の書き方の混入。
#
# **これは組版の代わりにはならない。**通っても組めるとは限らない。
# 組版できる環境では Makefile を使うこと。
#
# 使い方:
#   python3 check_tex.py                  # report/ の .tex 全部
#   python3 check_tex.py tacode_v2_en.tex

import argparse
import collections
import glob
import os
import re
import sys

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
DIRECTORY_FIGURE = os.path.join(DIRECTORY_SCRIPT, 'figure')
SUFFIX_FIGURE    = ('', '.pdf', '.png', '.jpg', '.eps')

# LaTeX に紛れ込みがちな Markdown の書き方。
# **箇条書きの検査は数式の外だけで行う**（align の中の行頭の `-` は引き算である）。
# コードスパンの検査は TeX の引用符 ``...'' を除くために、前後が backquote でない
# ものだけを見る
PATTERN_MARKDOWN = [
  (re.compile(r'\*\*'),            'Markdown bold (**) -- use \\textbf{}'),
  (re.compile(r'^#{1,6}\s', re.M), 'Markdown heading -- use \\section'),
  (re.compile(r'(?<!`)`(?!`)[^`\n]+(?<!`)`(?!`)'),
                                   'Markdown code span -- use \\texttt{}'),
]

PATTERN_MARKDOWN_TEXT = [
  (re.compile(r'^\s*[-*]\s+\S', re.M), 'Markdown list item -- use \\item'),
]

# 数式環境。ここは上の「本文だけの検査」から外す
ENVIRONMENT_MATH = ('align', 'align*', 'equation', 'equation*', 'gather', 'gather*',
                    'eqnarray', 'eqnarray*', 'array', 'bmatrix', 'pmatrix',
                    'lstlisting', 'verbatim')


def strip_comment(text):
  # 行末コメントを外す（\% は残す）
  out = []
  for line in text.split('\n'):
    index, escaped = None, False
    for position, character in enumerate(line):
      if escaped :
        escaped = False
        continue
      if character == '\\' :
        escaped = True
      elif character == '%' :
        index = position
        break
    out.append(line if index is None else line[:index])
  return '\n'.join(out)


def strip_math(text):
  # 数式環境と $...$ の中身を空白に置き換える（行番号は保つ）
  out = list(text)
  for name in ENVIRONMENT_MATH:
    pattern = re.compile(r'\\begin\{' + re.escape(name) + r'\}(.*?)\\end\{'
                         + re.escape(name) + r'\}', re.S)
    for match in pattern.finditer(text):
      for position in range(match.start(1), match.end(1)):
        if out[position] != '\n' :
          out[position] = ' '
  for match in re.finditer(r'(?<!\\)\$[^$\n]*(?<!\\)\$', ''.join(out)):
    for position in range(match.start(), match.end()):
      out[position] = ' '
  return ''.join(out)


def check(path):
  raw  = open(path, encoding='utf-8').read()
  text = strip_comment(raw)
  problem = []

  # 環境の対応
  depth = collections.Counter()
  for match in re.finditer(r'\\(begin|end)\{([a-zA-Z*]+)\}', text):
    depth[match.group(2)] += 1 if match.group(1) == 'begin' else -1
  for name, count in sorted(depth.items()):
    if count != 0 :
      problem.append('environment {:s} is unbalanced by {:+d}'.format(name, count))

  # 括弧と $。`\{` と `\}` は文字そのものなので数えない
  balance = len(re.findall(r'(?<!\\)\{', text)) - len(re.findall(r'(?<!\\)\}', text))
  if balance != 0 :
    problem.append('braces are unbalanced by {:+d}'.format(balance))
  if len(re.findall(r'(?<!\\)\$', text)) % 2 != 0 :
    problem.append('an odd number of $')

  # 参照
  label = set(re.findall(r'\\label\{([^}]+)\}', text))
  ref   = set(re.findall(r'\\(?:ref|eqref|pageref)\{([^}]+)\}', text))
  for name in sorted(ref - label):
    problem.append('\\ref{{{:s}}} has no \\label'.format(name))

  cite = {item.strip() for group in re.findall(r'\\cite\{([^}]+)\}', text)
          for item in group.split(',')}
  bib  = set(re.findall(r'\\bibitem\{([^}]+)\}', text))
  for name in sorted(cite - bib):
    problem.append('\\cite{{{:s}}} has no \\bibitem'.format(name))
  for name in sorted(bib - cite):
    problem.append('\\bibitem{{{:s}}} is never cited'.format(name))

  # 図
  for name in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text):
    if not any(os.path.exists(os.path.join(DIRECTORY_FIGURE, name + suffix))
               for suffix in SUFFIX_FIGURE) :
      problem.append('figure {:s} not found under figure/'.format(name))

  # Markdown の混入
  for pattern, message in PATTERN_MARKDOWN:
    for match in pattern.finditer(text):
      line = text[:match.start()].count('\n') + 1
      problem.append('line {:d}: {:s}'.format(line, message))

  outside = strip_math(text)
  for pattern, message in PATTERN_MARKDOWN_TEXT:
    for match in pattern.finditer(outside):
      line = outside[:match.start()].count('\n') + 1
      problem.append('line {:d} (outside maths): {:s}'.format(line, message))

  # \label が重複していないか
  duplicate = [name for name, count
               in collections.Counter(re.findall(r'\\label\{([^}]+)\}', text)).items()
               if count > 1]
  for name in sorted(duplicate):
    problem.append('\\label{{{:s}}} is defined more than once'.format(name))

  return problem


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('file', nargs='*', help='.tex files to check')
  args = parser.parse_args()

  file_list = args.file or sorted(glob.glob(os.path.join(DIRECTORY_SCRIPT, '*.tex')))
  flag_bad = False
  for path in file_list:
    problem = check(path)
    print('{:s}: {:s}'.format(os.path.basename(path),
                              'ok' if not problem else
                              '{:d} problem(s)'.format(len(problem))))
    for item in problem:
      print('  ' + item)
    flag_bad = flag_bad or bool(problem)

  print('')
  print('--This is not a substitute for typesetting: passing here does not mean it')
  print('  compiles. Use the Makefile where a TeX system is available.')
  sys.exit(1 if flag_bad else 0)


if __name__ == '__main__':
  main()

#set text(lang: "en", font: "KoPubWorldDotum_Pro", weight: "medium", size: 12pt)
#set page(margin: 2.0cm, paper: "a4")
#set page(numbering: "1")
#set math.equation(numbering: "(Eq. 1)", supplement: [Eq.])
#set figure(numbering: "1", supplement: [Fig.])

#show heading.where(level: 2): it => {
  pagebreak(weak: true)
  it
}

#show raw.where(block: true): set text(size: 8pt)
#let t(body) = highlight(fill: rgb("50C878"), body)
#let cf(body) = block(
  fill: rgb("F0F0F0"),
  inset: 10pt,
  above: 10pt,
  below: 10pt,
)[#text(weight: "light", fill: rgb("000000"))[c.f.] #h(1em) #body]

= Supplementary Text 1
Genomic sequence analysis using continuous wavelet transform

#outline(
  title: [Contents],
  target: heading.where(level: 2),
  indent: 1em,
)

== Signal transform of genomic sequences using complex numbers
=== Complex numbers
Complex numbers are represented as $z = p + i q$, while $i^2 = - 1$. $p$ is called as the real part, while $q$ is the imaginary part.
This represents a point $\( p \, q \)$ on the complex plane.
The size and conjugates of complex numbers are defined as follows


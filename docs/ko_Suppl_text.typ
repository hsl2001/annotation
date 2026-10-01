#set text(lang: "kr", font: "KoPubWorldDotum_Pro", weight: "medium", size: 12pt)
#set page(margin: 2.0cm, paper: "a4")
#set page(numbering: "1")
#show raw.where(block: true): set text(size: 8pt)
#let t(body) = highlight(fill: rgb("50C878"), body)

#t[supplementary text를 작성하기 위해 한국어로 미리 작성한 문서]

= Supplementary Text 1 (한국어 사전작성 버전)
- 연속 웨이블릿 변환을 이용한 유전체 서열 분석

#outline(
  title: [목차],
  target: heading.where(level: 2),
  indent: 1em,
)

== 1. 복소수를 통한 유전체의 신호 표현
=== 1.1. 복소수
- 복소수는 다음의 방식으로 표현 가능 $z = p + i q$, 이때 $i^2 = - 1$.
  - $p$는 실수부, $q$는 허수부.
  - 또는 평면상의 점 $\( p \, q \)$으로도 표현 가능 (복소평면).
  - 켤레복소수와 복소수 크기의 정의:
    $ overline(z) = p - i q \, #h(2em) \| z \| = sqrt(p^2 + q^2) \, #h(2em) z overline(z) = \| z \|^2 . $
  - 복소 함수의 적분은 실수부와 허수부를 따로 적분한 것.
    $ integral \( p \( t \) + i q \( t \) \) thin d t = integral p \( t \) thin d t + i integral q \( t \) thin d t . $


#t[여기 tmp]



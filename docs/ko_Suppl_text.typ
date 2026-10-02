#set text(lang: "kr", font: "KoPubWorldDotum_Pro", weight: "medium", size: 10pt)
#set page(margin: 2.0cm, paper: "a4")
#set page(numbering: "1")
#set math.equation(numbering: "(1)")

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

#t[supplementary text 한국어 버전]

= Supplementary Text 1
연속 웨이블릿 변환을 이용한 유전체 서열 분석

#outline(
  title: [목차],
  target: heading.where(level: 2),
  indent: 1em,
)

== 1. 복소수를 통한 유전체의 신호 표현
=== 1.1. 복소수
복소수는 $z = p + i q$, $i^2 = - 1$로 표현 가능하다. 이 때, $p$는 실수부, $q$는 허수부이다.
이는 복소평면상의 점 $\( p \, q \)$에 해당한다.
켤레복소수와 복소수의 크기는 다음과 같이 정의된다.
$
    overline(z) & = p - i q, \
        \| z \| & = sqrt(p^2 + q^2), \
  z overline(z) & = \| z \|^2 .
$ <complex-defs>
복소 함수의 적분은 실수부와 허수부를 따로 적분한 것이다.
$
  integral \( p \( t \) + i q \( t \) \) thin d t = integral p \( t \) thin d t + i integral q \( t \) thin d t .
$ <complex-integral-def>

=== 1.2. 회전과 복소수
오일러 공식에 의해,
$
  e^(i theta) = cos theta + i sin theta .
$ <euler-equation>
따라서 $e^(i omega t)$는 크기 1인 점이 각속도 $omega$로 회전하는 신호이다.
주파수 $f$는 '단위 길이당 회전 횟수'이고, 이는 $omega = 2 pi f$이다.
#cf[DNA에서 길이의 단위는  1bp로 한다.]

=== 1.3. DNA의 복소수 인코딩
기본값으로 아래의 인코딩을 사용한다.
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto, auto, auto),
    table.header([염기], [복소 값], [복소평면 좌표]),
    table.hline(),
    [A], [$1$], [$\( 1 \, 0 \)$],
    [C], [$i$], [$\( 0 \, 1 \)$],
    [G], [$- i$], [$\( 0 \, - 1 \)$],
    [T], [$- 1$], [$\( - 1 \, 0 \)$],
  )],
  kind: table,
)
이 때, 그 외 문자(`N` 등)는 0으로 처리한다.
예를 들어, `ACCTG`는 $1 \, i \, i \, - 1 \, - i$이다.
각 염기 값의 크기는 1이므로 신호 크기는 항상 1 이하라고 볼 수 있다.
#cf[위 인코딩은 생물학적으로 유일한 정의는 아니다. 이 인코딩에서 상보 염기의 값은 원래 값의 부호 반전이다. 즉, $z_(upright("complement")) = - z$이다. 다른 인코딩 방식을 선택하면 DNA 서열이 다른 신호로 변환되고, CWT 결과도 달라진다. ]

== 2. DNA 인코딩을 통한 연속시간 유전체 신호의 생성
=== 2.1. 연속시간 함수 표현의 필요성
연속 웨이블릿 변환은 짧은 진동 패턴(wavelet)이 타깃 신호와 얼마나 맞는지 비교하는 적분이다.
구체적으로, 신호와 켤레 wavelet의 복소 내적이다.
염기 번호 $n = 0 \, 1 \, dots.h \, L - 1$에 대응하는 값을 $x \[ n \]$이라고 표현할 수 있다.
이 때, 대괄호로 된 표현 $x \[ n \]$는 정수 위치에서만 정의된 배열이다.
같은 $x \[ 0 \] \, x \[ 1 \]$을 지나는 연속 위치의 함수는 여러 개이다.
두 점 사이를 직선으로 잇거나, 계단 모양으로도 만들 수 있다. 즉, 다양한 적분값이 나올 수 있다.
따라서 CWT 적분을 엄밀하게 정의하기 위해서는 불연속적인 점이 아니라 연속시간 함수를 정의해야 한다.
#cf[수학적으로 연속시간 함수(continuous-time function)와 연속함수(continuous function)는 다르다.
  연속시간 함수는 정의역이 실수 시간변수인 함수라는 뜻이고, 연속함수는 그래프가 끊기거나 뛰지 않고 이어지는 함수라는 뜻이다.
]

=== 2.2. 연속시간 함수로의 표현
연속시간 함수로 만들기 위하여 반올림 함수처럼 생각하고,
$
  x \( t \) = x \[ n \] quad upright("if ") n - 1 / 2 lt.eq t < n + 1 / 2 .
$ <dna-signal-def>
로 정의가 가능하다.
예를 들어,
$
  x \( t \) = cases(delim: "{", 1 \, & - 1 / 2 lt.eq t < 1 / 2 \,, i \, & 1 / 2 lt.eq t < 3 / 2 \,, 0 \, & upright("other positions").)
$
이다.
이 때, 서열 밖은 0으로 정의된다.
이 방식을 통하여 계산할 신호를 정확히 정의할 수 있다.

== 3. DNA 신호를 분석하기 위한 Wavelet
=== 3.1. 위치별 서열 신호 패턴의 분석
푸리에 분석:

#set text(lang: "kr", font: "KoPubWorldDotum_Pro", weight: "medium", size: 10pt)
#set page(margin: 2.0cm, paper: "a4")
#set page(numbering: "1")
#set math.equation(numbering: "(Eq. 1)", supplement: [Eq.])

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
$ <example-dna-signal>
이다.
이 때, 서열 밖은 0으로 정의된다.
이 방식을 통하여 계산할 신호를 정확히 정의할 수 있다.

== 3. DNA 신호를 분석하기 위한 Wavelet
=== 3.1. 위치별 서열 신호 패턴의 분석
푸리에 분석은 신호 전체에 $e^(i omega t)$가 얼마나 포함되어 있는지를 정량할 수 있다.
그러나, 어느 위치에서 해당 패턴이 나타났는지는 전체 비교로는 알기 어렵다.
Wavelet은 일정한 주기의 진동에 짧은 창 함수(window function)을 씌워 특정 위치의 패턴을 비교한다.
선택된 Morlet wavelet은 중심값이 크고 양쪽에서 빠르게 작아지는 가우스 함수를 채택한다.
$
  g \( u \) = e^(- u^2 \/ 2)
$ <gauss-window>
이 때 상수로서 진동 $e^(6 i u)$을 곱하면
$
  g \( u \) e^(6 i u) = e^(- u^2 \/ 2) \( cos 6 u + i sin 6 u \)
$ <gauss-window-omega-6>
형태가 된다.
이 wavelet은 실수부는 중심값이 큰 코사인 진동, 허수부는 중심값이 큰 사인 진동이다.
임의로 정한 Mother wavelet의 중심 각속도 $omega_0 &= 6$는 반드시 6일 필요는 없다.
같은 스케일에서 $omega_0$만 바꾸면 가우스 창의 폭은 그대로이고 그 촘촘함이 변한다.
#cf[$omega_0$가 크면 Morlet wavelet 코일이 더 촘촘해진다고 보면 된다.]
많은 신호분석 문헌에서 6이 경험적으로 좋은 시간-주파수 해상도를 보이는 것이 알려져 있다.
#cf[시간-주파수 해상도와 하이젠버그의 불확정성 원리: \
  #t[to be written...]
]
여기서 Mother wavelet이란, wavelet을 sliding하고 scaling하는 기준이 되는, 스케일 1·위치 0의 wavelet이다.
=== 3.2. 웨이블릿의 조건 1 - 영평균 조건
신호가 모든 위치에서  $x \( t \) = c$인 상수라고 하면 직관적으로 상수 패턴에는 진동이 없으므로 wavelet과 비교한 결과도 0이 되는 것이 적절하다.
그러나, 비교 패턴 자체의 적분이 0이 아니면 상수 신호에도 비교한 결과가 0이 나오지 않는다.
따라서 Mother wavelet(과 그 daugther wavelet들까지도) $psi \( u \)$는 다음 조건을 만족해야 한다.
$
  integral_(- oo)^oo psi \( u \) thin d u = 0 .
$ <zero-mean-condition>
Mother wavelet이 전체 구간에서의 정적분이 0이라는 것은 곧 평균값이 0이라는 뜻이다.
이를 #text(weight: "bold")[영평균 조건] 이라고 한다.

#cf[이외에도 유한 에너지 조건(finite energy condition)과 적합성 조건(admissibility condition)이 있다. 많은 Mother wavelet(충분히 매끄럽고 빠르게 감소하는 함수)이 적합성 조건을 만족하면 자동으로 영평균 조건을 만족한다.]

== 4. 영평균 조건을 만족하는 유전체 분석용 Morlet wavelet의 유도
=== 4.1. 가우스 적분
가우스 적분은 초월함수 적분법의 일종으로, 가우스 곡선의 $-oo$부터 $oo$까지의 정적분이다.
#cf[정규분포가 가우스 곡선의 일종이다.]
여기에서 사용할 두 적분은 다음과 같다.
$
  integral_(bb(R)) e^(- u^2 \/ 2) thin d u = sqrt(2 pi)
$ <gauss-integral-1>
$
  integral_(bb(R)) e^(- u^2 \/ 2) e^(i s u) thin d u = sqrt(2 pi) e^(- s^2 \/ 2) .
$ <gauss-integral-2>
여기서 $bb(R)$은 실수 전체를 뜻한다. 즉, 적분 구간은 $- oo$부터 $oo$까지이다.
첫 번째 식은 가우스 곡선 아래의 전체 면적을 의미한다.
두 번째 식은 가우스 함수에 각속도 $s$인 복소 진동을 곱한 적분이다.
이 때, $s = 0$이면 첫 번째 식과 같아지고, $s$가 커질수록 진동의 양과 음이 더 많이 상쇄되어 적분값이 작아진다.
=== 4.2. 가우스 적분의 유도 - 1
@gauss-integral-1 을 증명하기 위해, $J = integral_(bb(R)) e^(- u^2 \/ 2) thin d u$라고 놓을 수 있다.
같은 적분을 두 번 곱하며 두 적분의 변수를 $x$와 $y$로 구분하면,
$
  J^2 = integral_(bb(R)) integral_(bb(R)) e^(- \( x^2 + y^2 \) \/ 2) thin d x thin d y .
$ <gauss-integral-1-proof-1>
이는 평면 전체에서 높이 $e^(- \( x^2 + y^2 \) \/ 2)$를 면적에 곱해 더한 값이다.
높이가 원점으로부터의 거리에만 의존하므로 극좌표로 바꾸어 계산하면,
$x = r cos theta$, $y = r sin theta$로 놓으면 $x^2 + y^2 = r^2$이다.
이 때, $r$은 원점으로부터의 거리, $theta$는 각도이다.
반지름 폭 $d r$, 각도 폭 $d theta$인 작은 부채꼴의 면적은 호 길이 $r thin d theta$와 폭 $d r$의 곱인 $r thin d r thin d theta$이므로,
$
  J^2 = integral_0^(2 pi) integral_0^oo e^(- r^2 \/ 2) r thin d r thin d theta = 2 pi \[ - e^(- r^2 \/ 2) \]_0^oo = 2 pi .
$ <gauss-integral-1-proof-2>
$J > 0$이므로 $J = sqrt(2 pi)$이다.
#cf[이 적분은 수렴하고 적분항이 음수가 아님을 증명할 수 있으므로, 두 적분을 이중적분으로 묶을 수 있다.]
=== 4.3. 가우스 적분의 유도 - 2
이번에는 @gauss-integral-2 를 증명하기 위해,  각주파수 $s$에 따른 적분값을 $F \( s \)$라고 놓아보자.
$
  F \( s \) = integral_(bb(R)) e^(- u^2 \/ 2) e^(i s u) thin d u .
$ <gauss-integral-2-proof-1>
가우스 함수는 $\| u \|$를 곱해도 적분 가능할 정도로 빠르게 감소하는 것이 알려져 있다.
따라서, $s$에 대한 미분을 적분 안으로 옮길 수 있다.
$e^(i s u)$를 $s$로 미분하면 $i u e^(i s u)$이다.
#cf[이러한 교환이 가능한 이유:\
  #t[to be written...]
]
$
  F' \( s \) = i integral_(bb(R)) u e^(- u^2 \/ 2) e^(i s u) thin d u .
$ <gauss-integral-2-proof-2>
이 때, $g \( u \) = e^(- u^2 \/ 2)$이면 $g' \( u \) = - u g \( u \)$이므로 부분적분을 적용하면
$
  F' \( s \) & = - i integral_(bb(R)) g' \( u \) e^(i s u) thin d u \
             & = - i \[ g \( u \) e^(i s u) \]_(- oo)^oo + i integral_(bb(R)) g \( u \) \( i s \) e^(i s u) thin d u \
             & = - s F \( s \) .
$ <gauss-integral-2-proof-3>
이다.
경계항은 $\| e^(i s u) \| = 1$이고 $g \( u \) arrow.r 0$이므로 0이다.
따라서 $F' \( s \) = - s F \( s \)$이다.
곱의 미분법을 적용하면,
$
  frac(d, d s) \[ e^(s^2 \/ 2) F \( s \) \] = e^(s^2 \/ 2) \[ s F \( s \) + F' \( s \) \] = 0 .
$ <gauss-integral-2-proof-4>
따라서 $e^(s^2 \/ 2) F \( s \)$의 값은 상수이다.
$F \( 0 \) = sqrt(2 pi)$를 대입하면,
$
  F \( s \) = sqrt(2 pi) e^(- s^2 \/ 2) .
$ <gauss-integral-2-proof-4>
이다.
실수부는 코사인 적분 항등식이고, 허수부의 사인 적분은 0이다.
가우스 곡선은 우함수이고 사인함수는 기함수이므로 양쪽 기여가 상쇄된다.
$v = sqrt(2) u$로 치환하면 아래 식도 얻을 수 있다.
$
  integral_(bb(R)) e^(- u^2) e^(i s u) thin d u = sqrt(pi) e^(- s^2 \/ 4)
$ <gauss-integral-2-proof-5>
이다.
여기에 $s = 6$을 대입한 것이 $sqrt(pi) e^(- 9)$이다.
=== 4.4. 영평균 조건을 만족하는 보정항 계산
위의 wavelet, 즉 window function $times e^(i s u)$를 통해 다음 식을 계산할 수 있다.
$
  integral_(bb(R)) e^(- u^2 \/ 2) e^(i s u) thin d u = sqrt(2 pi) e^(- s^2 \/ 2)
$ <zero-mean-proof-1>
이 때, @zero-mean-proof-1 의 적분값은 0이 아니다.
따라서 $c$를 곱한 항을 빼서 $c e^(i s u)$ 적분값을 0으로 만들기 위해 시도해 볼 수 있다.
위의 각속도 $s = 6$을 대입하면,
$
  psi_0 \( u \) = e^(- u^2 \/ 2) \( e^(6 i u) - c \) .
$ <zero-mean-proof-2>
여기서 $psi_0$는 에너지 정규화 전의 wavelet임. 영평균 조건에 대입하면,
$
  0 = integral_(bb(R)) psi_0 \( u \) thin d u = sqrt(2 pi) e^(- 6^2 \/ 2) - c sqrt(2 pi) .
$ <zero-mean-proof-2>
양변을 $sqrt(2 pi)$로 나누면, $c = e^(- 18) .$
일반적인 각주파수 $omega_0$를 사용하면 보정값은 $e^(- omega_0^2 \/ 2)$이다.
$e^(- 18) approx 1.52 times 10^(- 8)$로 매우 작지만, 엄밀하게 평균 0을 만족시키기 위해 포함하였다.
#cf[현재 입력 DNA 신호의 평균을 빼는 것이 아니라 비교할 wavelet 자체를 보정하고 있다.]
== 5. 유한 에너지 조건을 만족하는 Morlet wavelet의 유도 (에너지 정규화)
=== 5.1. 에너지의 정의와 에너지 정규화의 목적
유한 에너지 조건은 wavelet 식의 제곱의 적분이 유한해야 한다는 조건이다.
에너지(제곱의 적분)이 무한하면 일반적인 유한한 신호와의 내적이 유한하다고 보장할 수 없다.
단, 여기서 정규화의 대상이 되는 것은 입력 DNA가 아니라 비교할 wavelet이다.
에너지 정규화라는 것은 직관적으로는 비교 패턴(wavelet)의 크기 기준을 고정하는 것이다.
$
  E = integral_(bb(R)) \| psi \( u \) \|^2 thin d u .
$ <energy-def>
어떤 함수의 전체 에너지는 @energy-def 와 같이 크기의 제곱을 전체 구간에서 적분한 값으로 정의할 수 있다.
이를 L2 에너지라고 부른다.
#cf[유클리드 거리, 즉 제곱합의 (제곱근)을 L2 거리(또는 노름; Norm)라고 부르기 때문이다. 즉, 복소평면에서 각 함숫값의 원점까지의 거리 제곱을 적분한 것이다. ]
이는 수식적으로는 관련이 있으나, 실제 물리적 에너지(운동에너지, 열에너지 등)나 생물학적 활성도를 뜻하는 것은 아니다.
wavelet에 상수 $N$을 곱하면 에너지는 $N^2$배가 된다.
$0 < E < oo$이면 $N = 1 \/ sqrt(E)$로 전체 에너지를 1로 맞출 수 있다.
이를 정규화(normalization)라고 부르며, 비교 패턴의 크기 기준을 고정하는 역할을 한다.
=== 5.2. 에너지 식의 전개
위에서 구한 영평균 조건을 위한 상수 $c = e^(- 18)$와 $psi \( u \) = N psi_0 \( u \)$라고 하자.
이 때, $N$은 양의 실수인 정규화 계수이다.
복소수 크기의 제곱은 자기 자신과 켤레 복소수의 곱이므로 식을 전개하면,
$
  \| e^(6 i u) - c \|^2 & = \( e^(6 i u) - c \) \( e^(- 6 i u) - c \) \
                        & = 1 + c^2 - c \( e^(6 i u) + e^(- 6 i u) \) \
                        & = 1 + c^2 - 2 c cos 6 u .
$ <energy-proof-1>
그러므로
$
  E = N^2 integral_(bb(R)) e^(- u^2) \( 1 + c^2 - 2 c cos 6 u \) thin d u .
$ <energy-proof-2>
@gauss-integral-2 에서  $v = sqrt(2) u$로 치환하고 실수부를 취하면,
$
  integral_(bb(R)) e^(- u^2) cos 6 u thin d u = sqrt(pi) e^(- 9) .
$ <energy-proof-2>
따라서
$
  E & = N^2 sqrt(pi) \( 1 + c^2 - 2 c e^(- 9) \) \
    & = N^2 sqrt(pi) \( 1 + e^(- 36) - 2 e^(- 27) \) .
$ <energy-proof-3>
$E = 1$을 요구하여 $N$을 구하면 최종적으로
$
  N = \[ sqrt(pi) \( 1 + e^(- 36) - 2 e^(- 27) \) \]^(- 1 \/ 2)
$ <energy-proof-4>
를 사용하면 전체 에너지를 1로 맞출 수 있다.
이 때, $e^(- 36)$은 $c^2$, $e^(- 27)$은 $c e^(- 9)$에서 나온 항이다.
보정항이 작으므로 $N$은 $pi^(- 1 \/ 4) approx 0.7511$에 매우 가깝다.
== 6. 적합성 조건을 만족하는 Morlet wavelet의 유도

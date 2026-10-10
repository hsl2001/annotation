#import "@preview/cetz:0.5.2"

#set text(lang: "kr", font: "KoPubWorldDotum_Pro", weight: "medium", size: 10pt)
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

#t[사전 작성된 supplementary text 한국어 버전]

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
== 5. 유한 에너지 조건과 적합성 조건을 만족하는 Morlet wavelet의 유도
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
=== 5.3. 적합성 조건을 만족하는 Morlet wavelet의 유도
푸리에 변환은
$
  hat(psi) \( xi \) = integral_(bb(R)) psi \( u \) e^(- i xi u) thin d u
$ <fourier-transform>
으로 정의된다.
가우스 적분을 각 항에 적용하면,
$
  hat(psi) \( xi \) = N sqrt(2 pi) [e^(- \( xi - 6 \)^2 \/ 2) - e^(- 18) e^(- xi^2 \/ 2)] .
$ <fourier-transform-wavelet>
이 때, $xi$는 Mother wavelet 좌표에서의 각주파수이다.
$xi = 0$에서는 두 항이 같으므로 $hat(psi) \( 0 \) = 0$이다. 평균 0을 주파수 영역에서 표현한 것이다.
이 식은 0 근처에서는 미분 가능하므로, $\| hat(psi) \( xi \) \| lt.eq C \| xi \|$인 상수 $C$를 선택할 수 있다.
멀리서는 가우스 함수와 같이 감소하고, 따라서 아래 적합성 적분이 유한하다.
$
  0 < integral_(bb(R)) frac(\| hat(psi) \( xi \) \|^2, \| xi \|) thin d xi < oo
$ <admissibility-condition>
이 성립한다.
0 근처에서 분자의 크기 제곱은 $\| xi \|^2$에 비례하는 상한을 가지므로 $\| xi \|$로 나눠도 적분항은 발산하지 않는다.
무한대에서는 가우스 항의 빠른 감소에 의해 적분이 유한하다.
== 6. Scaling과 sliding을 통한 CWT의 정의
=== 6.1. Wavelet의 scaling과 sliding
Mother wavelet을 $psi \( u \)$라 하고,
$
  psi_(a \, b) \( t \) = 1 / sqrt(a) psi (frac(t - b, a))
$ <daugther-wavelets>
스케일 $a > 0$와 위치 $b$를 사용해 @daugther-wavelets 와 같이 daughter wavelet들을 정의할 수 있다.
이 때, $t - b$는 wavelet의 중심을 위치 $b$로 이동 (sliding)하는 것이고, $\( t - b \) \/ a$는 패턴의 폭을 $a$배로 변경 (scaling) 하는 것이다.
$1 \/ sqrt(a)$는 폭을 바꾸더라도 에너지를 1로 유지하는 계수이다.
예를 들어, Mother wavelet의 위치 $u = 2$는 실제 신호에서 $t = b + 2 a$에 해당한다. $a = 4$이면 중심에서 8 bp 떨어진 위치, $a = 8$이면 16 bp 떨어진 위치이다.
$a$가 커질수록 더 넓은 구간과 더 긴 주기의 염기 패턴을 비교한다.
#cf[$a$는 창 함수의 폭을 바꾸는 배율이고, 실제 포함되는 염기의 개수와 정확히 같은 것응 아니다. 가우스 창 함수는 무한히 이어지며, 실제 컴퓨터의 계산에서는 $\| t - b \| lt.eq 8 a$인 범위만 남긴다. 예를 들어 $a = 4$이면 중심 양쪽 약 32 bp 이다. 위 절단 방식은 뒤에서 추가로 다룬다.]
#cf[Mother wavelet은 scaling 1,  sliding 0인 기준 wavelet으로 볼 수 있다.]
=== 6.2. scaling 과정에서 에너지를 1로 유지하기 위한 계수
폭을 늘리기만 하면 적분 구간도 늘어나 에너지가 $a$배로 커진다. 이를 상쇄하기 위해 진폭에 $1 \/ sqrt(a)$를 곱한다.
치환적분을 통해, $u = \( t - b \) \/ a$로 치환하면 $d t = a thin d u$이므로, 즉 @energy-scaling 에 의해
$
  integral_(bb(R)) \| psi_(a \, b) \( t \) \|^2 thin d t & = 1 / a integral_(bb(R)) lr(|psi (frac(t - b, a))|)^2 d t \
                                                         & = 1 / a integral_(bb(R)) \| psi \( u \) \|^2 a thin d u = 1 .
$ <energy-scaling>
스케일이나 슬라이딩 위치가 바뀌어도 비교 웨이블릿의 에너지는 1로 유지된다.
=== 6.3. 복소 내적을 통한 CWT의 계산
$
  W_x \( a \, b \) = integral_(bb(R)) x \( t \) overline(psi_(a \, b) \( t \)) thin d t = 1 / sqrt(a) integral_(bb(R)) x \( t \) overline(psi \( \( t - b \) \/ a \)) thin d t .
$ <cwt-def-inner-product>
@cwt-def-inner-product 가 현재 이 연구에서 사용하는 CWT의 정의이다.
이는 곧 위치 $b$ 주변의 신호와 스케일 $a$인 wavelet을 비교한 하나의 복소 계수를 얻는 것이다.
$x \( t \)$와 $overline(psi_(a \, b) \( t \))$를 곱하고 전체 위치에서 적분하는 방식이다.
이러한 비교를 복소 내적 (complex inner product)라고 부른다.
같은 위상(phase; 복소평면상 복소수의 각도)으로 맞는 성분은 보강되고, 맞지 않는 성분은 상쇄될 수 있다.
신호가 정확히 $x \( t \) = psi_(a \, b) \( t \)$이면 $W_x \( a \, b \) = integral \| psi_(a \, b) \( t \) \|^2 d t = 1$이다.
#cf[CWT가 웨이블릿과 유전체 신호의 '겹치는 면적'을 계산한다고 하는 것은 사실 직관적인 비유이다. 엄밀하게는, 두 그래프 사이의 기하학적 공통 면적을 구하는 것이 아니라 #text(weight: "bold")[신호와 켤레 wavelet의 곱을 적분]하는 것이다. 위치에 따라 이동시키며 비교하므로 상관(correlation) 형태이고, wavelet의 방향을 뒤집어 커널로 사용하면 합성곱(convolution)으로도 표현될 수 있다. 이 연구에서는 합성곱은 단순 복소내적이나 corrlation에 비해 계산하기 훨씬 간편하므로 convolution을 일차적으로 구한다.]
계수 $W = p + i q$를 해석하면, 크기 $\| W \| = sqrt(p^2 + q^2)$는 해당 패턴에 대한 wavelet 응답의 크기이다. 또, 위상(phase)는 복소평면상에서 계수의 각도로, 신호와 wavelet의 진동 위치 관계를 반영한다. 파워(power)는 $\| W \|^2 = p^2 + q^2$로 계산되는 크기를 제곱한 실수 값이다.
#cf[비교 패턴의 에너지가 1이라고 해서 모든 계수가 0과 1 사이인 것은 아니다. 입력 신호(유전체 신호)는 1로 정규화하지 않았으므로 DNA 신호에서도 $\| W \| > 1$이 나올 수 있다.] 따라서 CWT 계수를 상관계수나 확률처럼 0부터 1 사이의 유사도 점수로 해석하면 안 된다.
#cf[유한 에너지 조건은 계수의 크기에 영향을 준다. 코시-슈바르츠 부등식(Cauchy-Schwarz) 부등식을 적용하면, $ \| W_x \( a \, b \) \| lt.eq \| x \|_2 \| psi_(a \, b) \|_2 = \| x \|_2 . $, 그리고 $\| x \|_2 = sqrt(integral_(bb(R)) \| x \( t \) \|^2 thin d t)$이다. 길이 $L$인 DNA 신호는 $\| x \|_2^2 = sum_(n = 0)^(L - 1) \| x \[ n \] \|^2 lt.eq L$이므로 계수도 유한하다.]
=== 6.4. 알려진 신호의 CWT 응답 예시
결과를 직접 계산할 수 있는 신호를 통해 위 정의를 검산할 수 있다.
아래에서는 상수 신호, 정현파에서 검산 작업을 수행하였다.
==== 6.4.1. 상수 신호의 응답
무한히 이어진 상수 신호 $x \( t \) = c$에 대해,
$
  W_x \( a \, b \) = c sqrt(a) integral_(bb(R)) overline(psi \( u \)) thin d u = 0 .
$ <constant-signal-1>
영평균 조건에 의해 모든 스케일과 위치에서 계수가 0이다.
#cf[한 종류의 염기만 있는 유한 contig은 무한히 이어진 상수 신호와 다르다. contig 밖은 0으로 정의되므로 양 끝에 값의 변화가 생긴다.]
==== 6.4.2. 정현파 신호의 응답
계산이 용이한 연속시간 정현파 신호 $x \( t \) = e^(i omega t)$를 가정하고 $t = b + a u$를 대입하면,
$
  W_x \( a \, b \) = sqrt(a) thin e^(i omega b) integral_(bb(R)) e^(i a omega u) overline(psi \( u \)) thin d u .
$ <sinusoid-integral-1>
@sinusoid-integral-1 에 $overline(psi \( u \)) = N e^(- u^2 \/ 2) \( e^(- 6 i u) - e^(- 18) \)$를 대입하면, 가우스 함수에 곱하는 진동의 각주파수는 $a omega - 6$과 $a omega$가 된다.
가우스 적분을 두 번 사용하면,
$
  W_x \( a \, b \) = N sqrt(2 pi a) thin e^(i omega b) [e^(- \( a omega - 6 \)^2 \/ 2) - e^(- 18) e^(- \( a omega \)^2 \/ 2)] .
$ <sinusoid-integral-2>
여기서 $e^(i omega b)$는 위치에 따른 복소 위상을 나타내고, $e^(- \( a omega - 6 \)^2 \/ 2)$는 $a omega$가 6에 가까울 때 커진다.
따라서 고정 스케일에서 중심 주파수와 대표 주기는 아래와 같이 근사될 수 있다.
$
  omega_c approx 6 / a \, #h(2em) f_c approx frac(6, 2 pi a) \, #h(2em) T approx frac(2 pi a, 6) upright(" bp")
$ <sinusoid-integral-result>
$T = 1 \/ f_c$이므로 '스케일 $a$ = 주기 $a$ bp'는 아니다.
#cf[입력 주파수를 고정하고 스케일을 바꾸며 최대 응답을 찾으면 앞의 $sqrt(a)$도 변함. 따라서 $a omega = 6$이 정확한 최대 응답 위치라는 뜻은 아니다. 예를 들어, 주기 3 bp에 맞는 스케일은 약 $a = 9 \/ pi approx 2.86$이다. ]
#cf[일반적인 $omega_0$에서는 $omega_c approx omega_0 \/ a$이다. 동일한 $omega_c$에서 $omega_0$를 키우려면 $a$도 키워야 한다. 스케일 1에서는 $f_c approx 0.955$ cycles/bp로, 간격 1인 샘플의 Nyquist sampling 기준 $0.5$ cycles/bp보다 크다.]
==== 6.4.3 작은 스케일의 유효성
정수 샘플에서는 @integer-sample 을 만족한다.
$
  e^(i \( omega + 2 pi k \) n) = e^(i omega n) quad \( k upright("is integer") \)
$ <integer-sample>
서로 다른 연속 주파수가 같은 샘플을 만들 수 있고, 이를 엘리어싱 (aliasing)이라고 한다.
따라서 샘플만으로는 원래의 연속 주파수를 유일하게 알아낼 수 없다.
Nyquist 샘플링 정리에 의한 복원 보정이 필요하다.
이 연구에서는 원 신호를 복원하는 대신 아예 연속시간함수로 정의된 계단 형태의 함수를 선택했다.
==== 6.4.4 정현파 샘플의 계단 신호화
순수 연속 정현파(순수 사인파; sine wave)와 유전체 신호는 다르다.
예를 들어, `ACGT` 반복은 $x \[ n \] = e^(i pi n \/ 2)$이지만, 각 구간 안에서는 상수이다.
따라서 그 CWT에 연속 정현파에서의 공식을 그대로 적용할 수 없다.
유전체 신호는 다음 회전 패턴들의 합으로 표현할 수 있다.
$
  omega = pi \/ 2\, #h(2em) Omega_k = omega + 2 pi k, #h(2em) x \( t \) = sum_(k in bb(Z)) c_k e^(i Omega_k t) \, #h(2em) c_k = frac(sin \( Omega_k \/ 2 \), Omega_k \/ 2) .
$ <genome-signal-basis>
@genome-signal-basis 는 푸리에 급수로, 주기 함수를 여러 주파수의 진동으로 분해한 식이다.
이 때 $bb(Z)$는 정수 전체이고, $Omega_k$들은 샘플에서는 구분되지 않는 주파수이다.
$c_k$는 폭 1인 유전체 신호의 계단 구간 때문에 생가는 각 성분의 계수이다.
#cf[계단의 점프 지점에서는 급수가 양쪽 값의 평균으로 수렴하고, 따라서 그 점들의 값들은 CWT 적분에 영향을 주지 않는다.
  이를 증명하자면,
  $
    integral_(- 1 \/ 2)^(1 \/ 2) e^(- i Omega_k t) thin d t = frac(2 sin \( Omega_k \/ 2 \), Omega_k) = c_k .
  $ <genome-signal-proof-1>
  각 폭 1인 구간에서의 적분은 @genome-signal-proof-1 와 같이 표현된다.
]
길이 4인 한 주기에서 급수의 계수는, $1 / 4 integral_(- 1 \/ 2)^(7 \/ 2) x \( t \) e^(- i Omega_k t) d t$이다. 네 구간으로 나누고 $t = n + v$로 치환하면 $e^(i omega n) e^(- i Omega_k n) = 1$이므로
위 구간 적분이 네 번 더해지고 앞의 $1 \/ 4$와 상쇄된다. 그 외 길이 4의 푸리에 주파수에서는 네 복소수의 등비합이 0이 되어 계수도 0이다.
CWT는 신호에 대해 선형이고, 신호의 합을 변환한 결과는 각 변환 결과의 합이므로, 각 회전 패턴의 식에서
정수 위치 $b$에서는 $e^(i Omega_k b) = e^(i omega b)$이므로
$
  W_x \( a \, b \) = e^(i omega b) N sqrt(2 pi a) sum_(k in bb(Z)) c_k [e^(- \( a Omega_k - 6 \)^2 \/ 2) - e^(- 18) e^(- \( a Omega_k \)^2 \/ 2)] .
$ <genome-signal-proof-2>
가우스 항 때문에 먼 주파수의 기여가 빠르게 작아진다.
유한 합으로 근사할 수 있다.
#cf[위 식은 무한 신호의 결과이고, contig 등으로 유한한 신호에서는 경계 효과가 추가된다.]
== 7. 연속 적분과 커널 합
cwt는 신호와 켤레 wavelet의 곱을 적분한 것이다.
이를 일차적으로 빠르게 하기 위해, 각 염기의 기여를 계산해 더하는 식으로 바꾼다.
`AC`에서 수행해 보면 우선,  `A` 구간에서 $- 1 \/ 2 lt.eq t < 1 \/ 2$에서 값이 1이고
`C` 구간에서 $1 \/ 2 lt.eq t < 3 \/ 2$에서 값이 $i$이며 나머지 구간은 값이 0이다.
이 때, 스케일 $a$와 위치 $b$를 고정하고 아래와 같이 계산할 수 있다.
$
  W_(upright("AC")) \( a \, b \)
  = 1 / sqrt(a) \[ integral_(- 1 \/ 2)^(1 \/ 2) overline(psi \( \( t - b \) \/ a \)) d t
  + i integral_(1 \/ 2)^(3 \/ 2) overline(psi \( \( t - b \) \/ a \)) d t \] .
$ <genome-signal-example-1>
염기 $n$의 구간은 $\[ n - 1 \/ 2 \, n + 1 \/ 2 \)$이므로
$
  W_x \( a \, b \) = 1 / sqrt(a) sum_(n = 0)^(L - 1) x \[ n \] integral_(n - 1 \/ 2)^(n + 1 \/ 2) overline(psi \( \( t - b \) \/ a \)) thin d t .
$ <genome-signal-example-2>
@genome-signal-example-2 의 의미는 '각 염기 값' $times$ '그 구간의 wavelet 적분을 모두 더한 값' 이다.
치환적분을 통해 실제 위치 $t$를 Mother wavelet 좌표 $u$로 바꾸면,
$u = \( t - b \) \/ a$, 그리고 $\( n - b - 1 \/ 2 \) \/ a$와 $\( n - b + 1 \/ 2 \) \/ a$이다.
또, $d t = a d u$이기 때문에 앞의 계수는 $a \/ sqrt(a) = sqrt(a)$ 이다.
$
  W_x \( a \, b \) = sqrt(a) sum_(n = 0)^(L - 1) x \[ n \] integral_(\( n - b - 1 \/ 2 \) \/ a)^(\( n - b + 1 \/ 2 \) \/ a) overline(psi \( u \)) thin d u .
$ <genome-signal-integral>
위 @genome-signal-integral 와 같이 계산된다.
@genome-signal-integral 의 적분 구간은 $j = n - b$에만 의존한다.
따라서 같은 스케일, 같은 거리에서는 적분값도 같다.
이 때문에 거리별로 한 번씩 적분해 저장해 두고, 이 적분 재열을 커널(kernel)이라고 부른다.
$
  k_a \[ j \] = sqrt(a) integral_(\( j - 1 \/ 2 \) \/ a)^(\( j + 1 \/ 2 \) \/ a) overline(psi \( u \)) thin d u
$ <kernel-def>
커널이 준비되면 이제 각 위치에서 필요한 것은 염기 값 $times$ 거리별 커널 값의 합 뿐이다.
$
  W_x \( a \, b \) = sum_n x \[ n \] k_a \[ n - b \]
$ <cwt-def-by-kernel>
예를 들어, `AC`에서 중심을 한 칸 옮기면 각 염기의 상대 거리만 바뀐다.
#figure(
  table(
    columns: 4,
    table.header([비교 위치], [A의 거리와 기여], [C의 거리와 기여], [합계]),
    [$b = 0$], [$j = 0$: $k_a \[ 0 \]$], [$j = 1$: $i k_a \[ 1 \]$], [$k_a \[ 0 \] + i k_a \[ 1 \]$],
    [$b = 1$], [$j = - 1$: $k_a \[ - 1 \]$], [$j = 0$: $i k_a \[ 0 \]$], [$k_a \[ - 1 \] + i k_a \[ 0 \]$],
  ),
  kind: table,
)
커널을 이용하면 중심을 이동해도 다시 적분할 필요가 없기 때문에 같은 커널을 모든 위치에서 재활용한다.
== 8. 유한 커널과 수치 정밀도
커널을 실제로 컴퓨터가 계산할 수 있는 유한한 수준의 배열로 계산하기 위해 다음 두 가지를 정해야 한다.
첫 번째, 절단 범위는 wavelet의 꼬리를 어디에서 버릴지 정해야 한다.
두 번째, 수치 적분은 남긴 구간의 적분을 유한한 함숫값으로 근사하는 것이다.
이 때, 수치 적분은 구분구적법의 응용이다.
=== 8.1. Wavelet 꼬리의 절단
Wavelet은 중심에서 멀어질수록 함숫값이 작아지는, 꼬리가 긴 형태이다.
컴퓨터로는 무한한 배열을 저장할 수 없으므로 Mother wavelet의 좌표에서 $\| u \| > 8$인 부분을 버린다.
이를 염기 좌표로 표현하면 $\| t - b \| lt.eq 8 a$이다.
$
  r_a = ceil.l 8 a + 1 \/ 2 ceil.r \, #h(2em) M_a = 2 r_a + 1 .
$ <finite-kernel>
이 때, $r_a$는 배열의 반경, $M_a$는 왼쪽 $r_a$개 + 중심 1개 + 오른쪽 $r_a$개를 합한 길이이다.
$ceil.l v ceil.r$는 $v$ 이상인 가장 작은 정수이다.
예를 들어, $a = 4$이면 $r_a = ceil.l 32.5 ceil.r = 33$이고 $M_a = 67$이다.
이 때는 거리 $- 33$부터 $33$까지 저장한다.
바깥쪽 여분 구간은 실제 적분 범위와 겹치지 않이 0이 된다.
자연수 스케일에서는 같은 계산을 $r_a = 8 a + 1$, $M_a = 16 a + 3$으로 쓸 수 있다.
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto, auto, auto),
    table.header([스케일 $a$], [반경 $r_a$], [커널 길이 $M_a$]),
    table.hline(),
    [1], [9], [19],
    [4], [33], [67],
    [9], [73], [147],
    [64], [513], [1027],
  )],
  kind: table,
)
한 점의 값이 작다고 하여 wavelet 꼬리 전체의 적분도 작은 것은 아니다.
먼 구간까지 전부 더한 값의 기여를 확인해야 한다.
$\| x \( t \) \| lt.eq 1$을 이용하면 절단 전후의 계수 차이를 아래의 값으로 제한할 수 있다.
$
  E_(upright(t a i l)) lt.eq frac(N \( 1 + e^(- 18) \), 4) sqrt(a) thin e^(- 32)
$ <finite-kernel-limit-1>
이 식을 통해 절단 오차의 상한을 확인하면, $a = 4$에서 약 $4.76 times 10^(- 15)$, $a = 9$에서 약 $7.14 times 10^(- 15)$ 이하이다.
@finite-kernel-limit-1 를 증명하기 위해 아래와 같이 쓸 수 있다.
삼각부등식을 통해, $\| e^(6 i u) - e^(- 18) \| lt.eq 1 + e^(- 18)$이다.
$
  \| psi \( u \) \| lt.eq N \( 1 + e^(- 18) \) e^(- u^2 \/ 2) .
$ <finite-kernel-limit-2>
@finite-kernel-limit-2 로 둘 수 있고, $\| x \( t \) \| lt.eq 1$이며 양쪽 꼬리가 대칭이므로,
$
  E_(upright(t a i l)) & lt.eq sqrt(a) integral_(\| u \| > 8) \| psi \( u \) \| thin d u \
                       & lt.eq 2 N \( 1 + e^(- 18) \) sqrt(a) integral_8^oo e^(- u^2 \/ 2) thin d u .
$ <finite-kernel-limit-3>
$u gt.eq 8$에서 $1 lt.eq u \/ 8$을 이용하면

$
  integral_8^oo e^(- u^2 \/ 2) thin d u lt.eq 1 / 8 integral_8^oo u e^(- u^2 \/ 2) thin d u = e^(- 32) / 8
$ <finite-kernel-limit-4>
=== 8.2. 르장드르 다항식을 이용한 수치 적분
스케일 $a$와 거리 $j$에 해당하는 원래 적분 구간은
$\[\(j - 1 \/ 2\) \/ a \, \(j + 1 \/ 2\) \/ a\]$이다.
구현에서는 이 구간을 $\[ - 8 \, 8 \]$과 겹치는 부분으로 자른다.
겹치는 부분이 없으면 적분값은 0이다.
남은 구간의 길이를 $d > 0$이라 하면,
$m = ceil.l 2 d ceil.r$개의 작은 구간으로 나눈다.
따라서 각 구간의 길이 $d \/ m$은 최대 0.5이다.
예를 들어, $a = 4$, $j = 0$이면 구간은 $\[ - 0.125 \, 0.125 \]$이고 길이는 0.25이다.
이때 $m = ceil.l 0.5 ceil.r = 1$이므로 한 구간만 사용한다.
$a = 1$, $j = 0$이면 길이가 1이므로 두 구간으로 나눈다.
각 작은 구간 $\[ l \, h \]$에서 중심 $v = \( l + h \) \/ 2$와 반길이 $q = \( h - l \) \/ 2$를 둔다.
기준 구간 $\[ - 1 \, 1 \]$의 위치 $z$를 $u = v + q z$로 옮기면,
$
  integral_l^h f \( u \) thin d u
  approx q sum_(r = 1)^8 w_r \[ f \( v - q z_r \) + f \( v + q z_r \) \] .
$ <legendre-integral-1>
#figure(
  cetz.canvas(length: 1cm, {
    import cetz.draw: *
    set-style(stroke: (thickness: 0.8pt, cap: "round"))
    let roots = (
      0.0950125,
      0.2816036,
      0.4580168,
      0.6178762,
      0.7554044,
      0.8656312,
      0.9445750,
      0.9894009,
    )
    let sample-height(x) = 2 + 0.55 * calc.sin(0.8 * x) + 0.25 * calc.cos(1.7 * x)
    let chosen = roots.at(3)
    let left = 5 * (1 - chosen)
    let right = 5 * (1 + chosen)
    line((0, 0), (10.5, 0), mark: (end: "stealth"))
    line((0, 0), (0, 3.15), mark: (end: "stealth"))
    content((10.5, -0.08), $u$, anchor: "north-west")
    content((-0.12, 3.15), $f(u)$, anchor: "south-east")
    line(
      ..(
        for i in range(0, 101) {
          let x = i / 10
          ((x, sample-height(x)),)
        }
      ),
      stroke: rgb("2563EB") + 1.5pt,
    )
    for root in roots {
      let x-left = 5 * (1 - root)
      let x-right = 5 * (1 + root)
      let is-chosen = root == chosen
      let paint = if is-chosen { rgb("EA580C") } else { rgb("94A3B8") }
      let thickness = if is-chosen { 1.2pt } else { 0.55pt }
      line((x-left, 0), (x-left, sample-height(x-left)), stroke: paint + thickness)
      line((x-right, 0), (x-right, sample-height(x-right)), stroke: paint + thickness)
      circle((x-left, sample-height(x-left)), radius: if is-chosen { 0.075 } else { 0.045 }, fill: paint, stroke: none)
      circle(
        (x-right, sample-height(x-right)),
        radius: if is-chosen { 0.075 } else { 0.045 },
        fill: paint,
        stroke: none,
      )
    }
    line((0, -0.25), (0, -0.38))
    line((5, -0.25), (5, -0.38))
    line((10, -0.25), (10, -0.38))
    content((0, -0.4), $l$, anchor: "north")
    content((5, -0.4), $v$, anchor: "north")
    content((10, -0.4), $h$, anchor: "north")
    line((5, -0.85), (10, -0.85), mark: (end: "stealth"))
    line((5, -0.85), (0, -0.85), mark: (end: "stealth"))
    content((7.5, -0.72), $q$, anchor: "south")
    content((5, -1.05), [$q = (h-l)/2$], anchor: "north")
    for (x, label) in (
      (0, $-1$),
      (left, $-z_r$),
      (5, $0$),
      (right, $z_r$),
      (10, $1$),
    ) {
      line((x, -1.65), (x, 0), stroke: rgb("CBD5E1") + 0.55pt)
      line((x, -1.58), (x, -1.72), stroke: rgb("334155") + 0.75pt)
      content((x, -1.76), label, anchor: "north")
    }
    line((0, -1.65), (10.25, -1.65), mark: (end: "stealth"))
    content((10.3, -1.68), $z$, anchor: "west")
    content((5, -2.25), [기준 구간: $z in [-1,1]$], anchor: "north")
    content(
      (5, -2.75),
      [
        붉은색: $w_r [f(v-q z_r) + f(v+q z_r)]$ \
        근사된 전체 적분: $q sum_(r=1)^8$ (각 대칭쌍의 가중합)
      ],
      anchor: "north",
    )
  }),
  caption: [기준 구간의 르장드르 근 $plus.minus z_r$를 $u = v + q z$로 $[l,h]$에 옮긴 평가점이 주황색 점이다(@legendre-integral-1). 각 대칭쌍의 함수값을 더해 가중치 $w_r$를 곱하고, 치환에 따른 폭 계수 $q$를 마지막에 곱한다.],
)
#cf[
  르장드르 다항식 $P_n$은 구간 $\[ - 1 \, 1 \]$에서 낮은 차수의 다항식과 직교하는 다항식 계열이다.
  점화식으로 차수를 하나씩 높여 만들 수 있다.
  $
    P_0 \( z \) = 1 \, quad P_1 \( z \) = z \, quad
    P_(n + 1) \( z \) = frac(\( 2 n + 1 \) z P_n \( z \) - n P_(n - 1) \( z \), n + 1) .
  $ <legendre-integral-2>
  직교한다는 것은 $deg(q) < n$인 다항식 $q$에 대해(단, $deg$는 차수)
  $integral_(- 1)^1 P_n \( z \) q \( z \) thin d z = 0$임을 뜻한다.
  예를 들어 $P_2 \( z \) = \( 3 z^2 - 1 \) \/ 2$의 근은 $plus.minus 1 \/ sqrt(3)$이다.
  16점 구적법에서는 $P_16$의 서로 다른 실근 16개를 적분 위치로 사용하며, 이 근들은 모두 $\[ - 1 \, 1 \]$ 안에 있다.
  $P_16$은 짝수 차수 다항식이므로 근은 0을 기준으로 대칭이다.
  따라서 코드처럼 $z_r$인 양의 근 8개와 가중치만 저장하고, 각 $z_r$와 $- z_r$를 한 쌍으로 평가한다.
  $w_r$는 각 근에 대응하는 가중치이며, 대칭인 두 근에는 같은 가중치가 대응한다.
  이 위치와 가중치를 쓰는 16점 가우스-르장드르 구적법은 31차 이하 다항식을 정확히 적분한다.
  그 이유는 31차 이하의 다항식 $p$를 $P_16$으로 나누어
  $p = q P_16 + r$로 쓰면 $deg(q) lt.eq 15$이고 $deg(r) lt.eq 15$이기 때문이다.
  직교성에 따라 $q P_16$의 적분은 0이고, 모든 구적 위치에서도 $P_16$이 0이므로
  이 항의 구적합 역시 0이다.
  남은 $r$의 적분은 가중치를 정할 때 정확히 맞춘다.
]
구현을 위해 정규화 전 켤레 Morlet 함수
$f \( u \) = e^(- u^2 \/ 2) \( e^(- 6 i u) - e^(- 18) \)$를 각 위치에서 평가한다.
$N = 1 \/ sqrt(sqrt(pi) \(1 + e^(- 36) - 2 e^(- 27)\))$을 곱하고,
커널을 만드는 단계에서 다시 $sqrt(a)$를 곱한다.
== 9. 커널 합의 FFT를 통한 계산
위에서 커널을 준비하는 과정과 이후의 적분 과정을 다루었다.
그러나, 긴 유전체 서열의 모든 위치/스케일에서 단순히 이 연산을 반복하면 계산량이 폭증한다.
여기서는 합성곱(Convolution)이라는 연산을 FFT(Fast Fourier Transform)을 통해 가속하는 법을 다룬다.
=== 9.1. 합성곱(Convolution)
합성곱은 두 배열(함수 또는 벡터)에서 값의 곱을 만들어 특정 규칙으로 더하는 계산이다.
입력을 $X = \[ 1, 2, 3 \]$, 커널을 $H = \[ 10, 20 \]$, 배열 밖의 값을 0이라고 하자. 아래와 같이 계산한다.
#figure(
  table(
    columns: 3,
    table.header([출력 인덱스 $m$], [곱해서 더하는 항], [결과]),
    [0], [$1 times 10$], [10],
    [1], [$2 times 10 + 1 times 20$], [40],
    [2], [$3 times 10 + 2 times 20$], [70],
    [3], [$3 times 20$], [60],
  ),
  kind: table,
)
합성곱 결과는 $Y = \[ 10, 40, 70, 60 \]$이다.
각 출력에서 사용하는 입력 인덱스와 커널 인덱스의 합이 출력 인덱스 $m$이 됨을 주목하자.
일반 신호 위치 $b$에 대해 합성곱은, @convolution-1 와 같이 작성된다.
$
  \( x * h \) \[ b \] = sum_n x \[ n \] h \[ b - n \]
$ <convolution-1>
앞서 원래 적분식 @genome-signal-integral 에서 각 염기 $n$의 적분 구간은 $n-b$에만 의존한다고 보았으므로, CWT 커널 인덱스 $n-b$는 염기 위치와 비교 위치 사이의 상대 거리를 뜻한다.
합성곱의 커널 인덱스는 $b - n$, CWT의 커널 인덱스는 $n - b$이다. 부호가 반대이므로, 아래와 같이 커널의 좌우를 뒤집는다.
$
  h_a \[ j \] = k_a \[ - j \]
$ <convolution-2>
스케일 $a$를 고정하고, CWT의 위치 $b$에서 입력값 $x \[ n \]$에 곱해지는 커널 값은 $k_a \[ n - b \]$이다.
반면 합성곱은 출력 위치 $b$에서 입력 $x \[ n \]$에 $h_a \[ b - n \]$을 곱한다.
두 식의 차이는 커널 인덱스의 부호뿐이다. @convolution-2 에 따라 합성곱에 사용할 커널을 $k_a$의 좌우 반전으로 정의하면,
$
  \( x * h_a \) \[ b \] & = sum_n x \[ n \] h_a \[ b - n \] \
                        & = sum_n x \[ n \] k_a \[ - \( b - n \) \] \
                        & = sum_n x \[ n \] k_a \[ n - b \] \
                        & = W_x \( a \, b \) .
$
즉, 각 입력 $x \[ n \]$과 그 입력 위치에서 $b$까지의 상대 거리 $n-b$에 해당하는 CWT 커널 값을 곱해 더하는 것과,
좌우를 뒤집은 커널 $h_a$로 합성곱을 계산하는 것은 항별로 완전히 같다.
따라서 스케일마다 $h_a \[ j \] = k_a \[ -j \]$를 한 번 준비하면, 합성곱의 출력 위치 $b$가 CWT의 위치 $b$에 대응하여 모든 위치의 CWT 계수를 구할 수 있다.
입력 배열 범위 밖의 $x \[ n \]$은 0으로 두므로, 경계에서도 두 계산의 합 범위가 일치한다.
#cf[음수 인덱스의 저장\
  거리 $j$는 음수일 수 있으나 저장 배열은 0부터 시작한다. 
  예를 들어, 반경 $r = 1$인 커널은 아래와 같이 저장된다. 
  #figure(
    table(
      columns: 4,
      table.header([배열 인덱스 $q$], [0], [1], [2]),
      [거리 $j$], [$- 1$], [0], [1],
      [저장 배열 $K$], [$k_a \[ - 1 \]$], [$k_a \[ 0 \]$], [$k_a \[ 1 \]$],
      [뒤집은 배열 $H$], [$k_a \[ 1 \]$], [$k_a \[ 0 \]$], [$k_a \[ - 1 \]$],
    ),
    kind: table,
  )
  반경 $r$에서는 거리 $- r, dots.h, r$를 인덱스 $0, dots.h, 2 r$에 저장한다. 즉, $K \[ q \] = k_a \[ q - r \]$이고 뒤집힌 배열은 $H \[ q \] = K \[ 2 r - q \] = k_a \[ r - q \]$이다. 중심이 배열의 $r$번에 있으므로 합성곱의 결과를 분석할 때도 그만큼의 이동을 반영해야 한다. 입력 배열이 contig의 위치 $s$부터 시작한다면, 원하는 위치 $b$의 계수는 입력 안의 위치 $b -s + r$ ($r$은 커널 중심의 배열 위치)에서 읽을 수 있다. 따라서 $m = b - s + r$와 같이 쓸 수 있다. 예를 들어 입력을 위치 100부터 읽었고, 위치 105의 계수가 필요하며, 커널 반경이 9이면 $m = 105 - 100 + 9 = 14$이다. 
  위 식 $m = b - s + r$을 증명하자면, 입력 배열의 상대 위치 $p$에 대해 $X \[ p \] = x \[ s + p \]$로 놓고 $Y \[ m \] = sum_p X \[ p \] H \[ m - p \] = sum_p x \[ s + p \] k_a \[ r - m + p \].$라 할 수 있다. 원하는 커널 인덱스가 $s + p - b$이므로 $r - m + p = s + p - b$를 풀면 $m = b - s + r$이다. 
]
=== 9.2. DFT를 통한 신호 패턴의 분해
DFT(Discrete Fourier Transform; 이산 푸리에 변환)는 위치별 값을 주파수별 값으로 바꾼다. 
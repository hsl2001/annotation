#set text(lang: "kr", font: "KoPubWorldDotum_Pro", weight: "medium", size: 12pt)
#set page(margin: 1.5cm, paper: "a4")
#set page(numbering: "1")
#show raw.where(block: true): set text(size: 8pt)
#let t(body) = highlight(fill: rgb("50C878"), body)

= 연속 웨이블릿 변환을 이용한 유전체 분석 방법
- DNA에서부터 유전체의 연속 웨이블릿 변환(continuous wavelet transform, CWT)으로 가는 과정.
- 최대한 모든 유도과정 및 계산 포함.
- #text(weight: "light")[생물학적 과정이 아닌 기술적인 부분만 다룸.]

== 내용
+ DNA를 복소수 좌표계에 매핑하고, 그 복소수 숫자들을 통해 유전체 신호로 만드는 방법.
+ 평균 0, 에너지 1인 Morlet mother wavelet을 만드는 방법.
+ Wavelet을 sliding하고 scaling하여 CWT를 정의하는 방법.
+ CWT 적분을 염기별 커널 합으로 바꾸고, 커널을 수치적으로 계산하는 방법.
+ 커널 합을 고속 푸리에 변환(FFT)으로 계산하는 방법.
+ 출력 행렬의 위치,스케일,파워를 해석하는 방법.

== 1. 복소수를 통해 유전체 신호를 표현
=== 1.1 복소수는 두 실수의 묶음
- 복소수 $z = p + i q$에서 $i^2 = - 1$.
- $p$는 실수부, $q$는 허수부.
- 이를 평면의 점 $\( p \, q \)$로 생각하면 됨(복소평면).
켤레복소수와 복소수 크기의 정의
$ overline(z) = p - i q \, #h(2em) \| z \| = sqrt(p^2 + q^2) \, #h(2em) z overline(z) = \| z \|^2 . $
복소 함수의 적분은 실수부와 허수부를 따로 적분한 것.
$ integral \( p \( t \) + i q \( t \) \) thin d t = integral p \( t \) thin d t + i integral q \( t \) thin d t . $
=== 1.2 회전을 복소수 지수로 표현 가능
오일러 공식:\
$ e^(i theta) = cos theta + i sin theta . $
따라서 $e^(i omega t)$는 크기 1인 점이 각속도 $omega$로 회전하는 신호.\
주파수 $f$를 “단위 길이당 회전 횟수”로 표현하면 $omega = 2 pi f$. \
c.f.) DNA에서 길이 단위는 1 bp로 둠.
=== 1.3 DNA의 숫자 표현은 선택 사항입니다
현재 구현에서는 다음 encoding을 사용.\
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
이 때, 그 외 문자(`N` 등)는 0으로 처리. \
e.g.) `ACCTG`는 $1 \, i \, i \, - 1 \, - i$. \
각 염기 값의 크기는 1이므로 신호 크기는 항상 1 이하임.
c.f.)\
위 encoding은 생물학적으로 유일한 정의는 아님. \
상보관계를 켤례관계(conjugate)로 표현하기는 하였음. => 그러나 수학적/생물학적 근거 없음. \
다른 encoding을 선택하면 DNA 신호가 다른 숫자 신호로 되고 CWT 결과도 달라짐.

== 2. 인코딩으로부터 연속인 신호 생성
=== 2.1 샘플링만으로는 적분할 함수를 알 수 없음
#text(
  weight: "light",
)[아래 내용은 유전체의 CWT를 정확히 정의하기 위해 2026년 하반기에 수학적 엄밀성을 위해 새로 도입한 개념임.]\
- #text(
    weight: "bold",
  )[중요한 사항) 연속 웨이블릿 변환을 요약하자면, 짧은 펄스 패턴(wavelet)가 타깃 신호와 겹치는 면적(엄밀하게는 합성곱; convolution)을 계산하는 것임. 즉, 적분임.]\
- 염기 번호 $n = 0 \, 1 \, dots.h \, L - 1$에 대응하는 값을 $x \[ n \]$이라고 하자.
- 이 때 대괄호로 된 표현 $x \[ n \]$는 정수 위치에서만 정의된 배열.
- 반면 $x \( t \)$의 소괄호는 실수 위치 $t$에서도 정의된 함수.
- 같은 $x \[ 0 \] \, x \[ 1 \]$을 지나는 연속 위치의 함수는 여러 개임.
  - 두 점 사이를 직선으로 잇거나, 계단 모양으로 만들 수도 있음.
  - 다양한 적분값이 나올 수 있음.
  - 따라서 CWT 적분을 엄밀하게 정의하기 위해서는 불연속적인 점이 아니라 연속시간 함수를 정의해야 함.
    - c.f.) 수학적으로 연속시간 함수(continuous-time function)와 연속함수(continuous function)는 다름.
    - 연속시간 함수는 정의역이 실수 시간변수인 함수라는 뜻이고, 연속함수는 그래프가 끊기거나 뛰지 않고 이어지는 함수.

=== 2.2 각 염기는 폭 1인 구간을 차지

$ x \( t \) = x \[ n \] quad upright("if ") n - 1 / 2 lt.eq t < n + 1 / 2 . $
반올림 함수라고 생각할 수 있음.\
예를 들어 `AC`라면
$
  x \( t \) = cases(delim: "{", 1 \, & - 1 / 2 lt.eq t < 1 / 2 \,, i \, & 1 / 2 lt.eq t < 3 / 2 \,, 0 \, & upright("other positions").)
$
이 때 서열 밖은 0으로 정의함.\
이 방식으로 계산할 신호를 정확하게 정의할 수 있음.\
== 3. Wavelet 만들기
=== 3.1 위치별 서열신호 패턴의 분석
푸리에 분석: \
신호 전체에 $e^(i omega t)$가 얼마나 포함되어있는지 비교 \
=> 그러나 어느 위치에서 그 패턴이 나타났는지는 전체 비교값 하나로 알기는 어려움.\
웨이블릿은 진동에 짧은 창 함수(window function)을 씌워 특정 위치 주변의 패턴을 비교함.\
여기에서 선택한 Morlet wavelet은 가운데가 크고 양쪽에서 빠르게 작아지는 가우스 함수를 창 함수로 선택함.
$ g \( u \) = e^(- u^2 \/ 2) . $
$g \( 0 \) = 1$이고 $g \( 2 \) = e^(- 2)$, $g \( 4 \) = e^(- 8)$임.
진동 $e^(6 i u)$를 곱하면
$ g \( u \) e^(6 i u) = e^(- u^2 \/ 2) \( cos 6 u + i sin 6 u \) $
실수부는 가운데가 큰 코사인 진동, 허수부는 가운데가 큰 사인 진동.\
- 이 때, 숫자 6은 구현에서 선택한 Mother wavelet의 각주파수 $omega_0$. 모든 CWT가 반드시 6을 사용해야 하는 것은 아님. \
- Mother wavelet의 각주파수 $omega_0$가 커지면 진동이 많아져 주파수 성분을 상대적으로 더 세밀하게 구분하고, 시간상 위치를 구분하는 능력은 상대적으로 낮아짐. 작아지면 그 반대.
- 다수의 신호 분석 문헌에서 6이 경험적으로 좋은 시간-주파수 해상도를 보이는 것이 알려져 있음.
c.f.) Mother wavelet이란: \
Wavelet을 sliding하고 scaling하는 기준점이 되는, 원점의 가장 짧은 wavelet.
=== 3.2 상수 신호에 의한 효과를 제거하기 위한 웨이블릿의 조건
부정적분을 하면 적분상수가 나옴. \
이 상수를 0으로 만들어야 하는 이유.
즉, $[-infinity \, infinity]$에서 적분했을 때 값이 0이어야 하는 이유?\
신호가 모든 위치에서 $x \( t \) = c$인 상수라고 하자.\
직관적으로, 상수 신호에는 진동이 없으므로 wavelet과 비교한 결과도 0이 되는 것이 적절함.\
그러나 비교 패턴의 적분이 0이 아니면 상수 성분에도 반응하게 됨. (비교한 결과가 0이 아님)\
따라서 Mother wavelet $psi \( u \)$에는 다음 조건을 부여함.\
$ integral_(- oo)^oo psi \( u \) thin d u = 0 . $
문장으로 풀어 쓰면, 웨이블릿의 평균값이 0이라는 뜻.\
곧, 실수부와 허수부의 적분이 각각 0이어야 한다는 뜻임.\
이는 wavelet의 조건 중 #text(weight: "bold")[평균 0 조건.]\
무한 구간 적분은 $lim_(R arrow.r oo) integral_(- R)^R$로 생각할 수 있음.
가우스 함수가 빠르게 감소하므로 여기서 사용하는 적분은 수렴함.
c.f.) 평균 0만으로 모든 함수가 wavelet이 되는 것은 아님. \
총 3가지 조건을 만족해야 하고, \
평균 0 조건(zero-mean condition) 이외에도 적합성 조건(admissibility condigion)과 유한 에너지 조건(finite energy condition)도 만족해야 함.
== 4. 평균이 0인 Morlet wavelet 유도하기
=== 4.1 가우스 적분
계산에 사용할 두 적분은 다음과 같음.
$
  integral_(bb(R)) e^(- u^2 \/ 2) thin d u = sqrt(2 pi) \, #h(2em) integral_(bb(R)) e^(- u^2 \/ 2) e^(i s u) thin d u = sqrt(2 pi) e^(- s^2 \/ 2) .
$
- $bb(R)$은 실수 전체를 뜻함. 즉, 적분 구간은 $- oo$부터 $oo$까지임.
- 첫 번째 식은 가우스 곡선 아래의 전체 면적.
- c.f.) 정규분포가 가우스 곡선의 일종임.
- 두 번째 식은 가우스 함수에 각주파수 $s$인 복소 진동을 곱한 적분.
  - $s = 0$이면 첫 번째 식과 같아짐.
  - $s$가 커질수록 진동의 양과 음이 더 많이 상쇄되어 적분값이 작아짐.
두 식의 유도는 뒤에서 다룸. 일종의 공식처럼 외우는 값임.\
지금은 이 결과를 이용해 보정항과 정규화 계수를 계산함.
=== 4.2 평균 0 조건으로 보정항 계산
위의 가우스 진동 패턴(창 함수 $times e^(i s u)$),\
즉, $integral_(bb(R)) e^(- u^2 \/ 2) e^(i s u) thin d u = sqrt(2 pi) e^(- s^2 \/ 2)$ 은 적분값이 0이 아님.\
따라서 가우스 창에 상수 $c$를 곱한 항을 빼서 적분값을 0으로 만들고자 함.\
위에서 제안한 각속도 $s = 6$을 대입하면
$ psi_0 \( u \) = e^(- u^2 \/ 2) \( e^(6 i u) - c \) . $
여기서 $psi_0$는 에너지 정규화 전의 wavelet임. 평균 0 조건에 대입하면
$ 0 = integral_(bb(R)) psi_0 \( u \) thin d u = sqrt(2 pi) e^(- 6^2 \/ 2) - c sqrt(2 pi) . $
양변을 $sqrt(2 pi)$로 나누면
$ c = e^(- 18) $
즉, `exp(-18.0)`는 #text(weight: "bold")[평균 0 조건에서 나온 값]임.\
일반적인 각주파수 $omega_0$를 사용하면 보정값은 $e^(- omega_0^2 \/ 2)$임.\
$e^(- 18) approx 1.52 times 10^(- 8)$로 매우 작지만, 엄밀하게 평균 0을 만족시키기 위해 포함함.\
주의할 점은 입력 DNA 신호의 평균을 빼는 것이 아니라 #text(weight: "bold")[비교할 wavelet 자체를 보정]하는 것임.
== 5. 에너지를 1로 맞추기
=== 5.1 에너지의 정의와 정규화 목적
앞서 서술한 유한 에너지 조건에 해당하는 부분.\
함수의 에너지는 크기의 제곱을 전체 구간에서 적분한 값으로 정의함.

$ E = integral_(bb(R)) \| psi \( u \) \|^2 thin d u . $

이를 연속 L2 에너지라고 부름. 실제 물리적 에너지나 생물학적 활성도를 뜻하는 것은 아님.
패턴의 전체 크기를 비교하기 위한 수학적 척도임.

wavelet에 상수 $N$을 곱하면 에너지는 $N^2$배가 됨.
따라서 $N$을 적절하게 선택해 에너지를 1로 맞출 수 있음.
이를 정규화(normalization)라고 부르며, 비교 패턴의 크기 기준을 고정하는 역할을 함.

=== 5.2 크기의 제곱 전개
<크기의-제곱을-전개합니다>
$c = e^(- 18)$로 두고, $psi \( u \) = N psi_0 \( u \)$로 놓음.
$N$은 양의 실수인 정규화 계수임.
복소수 크기의 제곱은 자기 자신과 켤레복소수의 곱이므로

$
  \| e^(6 i u) - c \|^2 & = \( e^(6 i u) - c \) \( e^(- 6 i u) - c \) \
                        & = 1 + c^2 - c \( e^(6 i u) + e^(- 6 i u) \) \
                        & = 1 + c^2 - 2 c cos 6 u .
$

그러므로

$ E = N^2 integral_(bb(R)) e^(- u^2) \( 1 + c^2 - 2 c cos 6 u \) thin d u . $

4절의 가우스 적분에서 변수 크기를 바꾸면

$
  integral_(bb(R)) e^(- u^2) thin d u = sqrt(pi) \, #h(2em) integral_(bb(R)) e^(- u^2) cos 6 u thin d u = sqrt(pi) e^(- 9) .
$

따라서

$
  E & = N^2 sqrt(pi) \( 1 + c^2 - 2 c e^(- 9) \) \
    & = N^2 sqrt(pi) \( 1 + e^(- 36) - 2 e^(- 27) \) .
$

$E = 1$을 요구하여 $N$을 구하면 최종적으로

$ #box(stroke: black, inset: 3pt, [$ N = \[ sqrt(pi) \( 1 + e^(- 36) - 2 e^(- 27) \) \]^(- 1 \/ 2) $]) $

임. $e^(- 36)$은 $c^2$, $e^(- 27)$은 $c e^(- 9)$에서 나온 항임.
보정항이 작으므로 $N$은 $pi^(- 1 \/ 4) approx 0.7511$에 매우 가까움.
그러나 구현에서는 위의 보정항을 포함한 식을 사용함.

c.f.) 이 정규화는 무한 구간에 정의된 #text(weight: "bold")[Mother wavelet의 연속 에너지]를 1로 맞추는 것임. \
나중에 구간 적분으로 만든 커널 배열의 제곱합은 반드시 1일 필요가 없음.
그 배열을 다시 정규화하면 원래 정의와 다른 변환을 계산하게 되므로 추가 정규화는 하지 않음.

== 6. 스케일과 위치로 CWT 정의하기
<스케일과-위치로-cwt-정의하기>
=== 6.1 Wavelet의 scaling과 sliding
<패턴을-늘리고-이동합니다>
Mother wavelet을 $psi \( u \)$라 하고, 스케일 $a > 0$와 위치 $b$를 사용해

$ psi_(a \, b) \( t \) = 1 / sqrt(a) psi #h(-1em) (frac(t - b, a)) $

로 정의함.

- $t - b$: 패턴의 중심을 위치 $b$로 이동(sliding).
- $\( t - b \) \/ a$: 패턴의 폭을 $a$배로 변경(scaling).
  - $a > 1$이면 늘리고, $0 < a < 1$이면 줄이는 것임.
- $1 \/ sqrt(a)$: 폭을 바꾸더라도 에너지를 1로 유지하는 계수.

e.g.) Mother wavelet의 위치 $u = 2$는 실제 신호에서 $t = b + 2 a$에 해당함. \
$a = 4$이면 중심에서 8 bp 떨어진 위치, $a = 8$이면 16 bp 떨어진 위치임.
$a$가 커질수록 더 넓은 구간과 더 긴 주기의 패턴을 비교함.
스케일은 주기나 커널 배열 길이와 같은 값은 아님. 주기와의 관계는 12.2절에서 계산함.

c.f.) Mother wavelet은 스케일 1, 위치 0의 기준 wavelet임. \
수학적으로 $a < 1$도 가능하므로 반드시 가장 짧은 wavelet이라는 뜻은 아님.

=== 6.2 Scaling 후의 에너지 유지
<왜-앞에-1sqrt-a를-붙이나요>
폭을 늘리기만 하면 적분 구간도 늘어나 에너지가 $a$배로 커짐.
이를 상쇄하기 위해 진폭에 $1 \/ sqrt(a)$를 곱함.
$u = \( t - b \) \/ a$로 치환하면 $d t = a thin d u$이므로

$
  integral_(bb(R)) \| psi_(a \, b) \( t \) \|^2 thin d t & = 1 / a integral_(bb(R)) lr(|psi #h(-1em) (frac(t - b, a))|)^2 d t\
  & = 1 / a integral_(bb(R)) \| psi \( u \) \|^2 a thin d u = 1 .
$

즉, 스케일이나 위치가 바뀌어도 비교 패턴의 에너지는 1로 유지됨.
$1 \/ a$를 사용하면 에너지 대신 다른 크기 기준을 따르는 변환이 됨.
따라서 다른 라이브러리의 결과와 비교할 때는 정규화 규칙도 확인해야 함.

=== 6.3 복소 내적으로 CWT 계산
<신호와-패턴의-겹침을-계산합니다>
두 복소 함수를 비교하기 위해 다음 적분을 사용함.

$
  #box(stroke: black, inset: 3pt, [$ W_x \( a \, b \) = integral_(bb(R)) x \( t \) overline(psi_(a \, b) \( t \)) thin d t $]) = 1 / sqrt(a) integral_(bb(R)) x \( t \) overline(psi \( \( t - b \) \/ a \)) thin d t .
$

이것이 현재 구현에서 사용하는 CWT 정의임.
위치 $b$ 주변의 신호와 스케일 $a$인 wavelet을 비교한 하나의 복소 계수를 얻는 것임.

- $x \( t \)$와 $overline(psi_(a \, b) \( t \))$를 곱하고 전체 위치에서 적분함.
- 이런 비교를 복소 내적(inner product)이라고 부름.
- 같은 위상으로 맞는 성분은 더해지고, 맞지 않는 성분은 상쇄될 수 있음.

신호가 정확히 $x \( t \) = psi_(a \, b) \( t \)$이면
$W_x \( a \, b \) = integral \| psi_(a \, b) \( t \) \|^2 d t = 1$임.
켤레를 사용하는 이유가 이 관계에서 드러남.
예를 들어 1.1절처럼 $i$와 $i$를 그대로 곱하면 $- 1$이지만,
$i$와 그 켤레 $- i$를 곱하면 1이 됨.

c.f.) “겹치는 면적”은 직관적인 비유임. \
두 그래프 사이의 기하학적 공통 면적을 구하는 것이 아니라, #text(weight: "bold")[신호와 켤레 wavelet의 곱을 적분]하는 것임.
위치에 따라 이동시키며 비교하므로 상관(correlation) 형태이고,
wavelet의 방향을 뒤집어 커널로 사용하면 합성곱(convolution)으로도 표현할 수 있음. 이 관계는 9.1절에서 다룸.

계수 $W = p + i q$의 해석:
- 크기 $\| W \| = sqrt(p^2 + q^2)$: 해당 패턴에 대한 응답의 크기.
- 위상: 복소평면에서 계수의 각도. 신호와 wavelet의 진동 위치 관계를 반영함.
- 파워(power) $\| W \|^2 = p^2 + q^2$: 크기를 제곱한 실수 값.

비교 패턴의 에너지가 1이라고 해서 모든 계수가 0부터 1 사이인 것은 아님.
입력 신호의 에너지는 1로 정규화하지 않았으므로 DNA 신호에서도 $\| W \| > 1$이 나올 수 있음.
따라서 CWT 계수를 확률이나 0부터 1 사이의 유사도 점수로 해석하면 안 됨.

=== 6.4 연속 CWT와 출력 배열의 구분
<연속-cwt와-출력-배열은-구분합니다>
수학적 정의에서는 $a > 0$와 $b$가 연속적으로 변할 수 있음.
그러나 컴퓨터에서 모든 실수 스케일과 위치의 값을 저장할 수는 없음.
현재 프로그램은 사용자가 선택한 자연수 스케일과 정수 위치에서만 CWT를 계산함.

즉, 출력은 #text(weight: "bold")[연속 CWT를 특정 스케일과 위치에서 샘플링한 배열]임.
계수를 배열에 저장한다고 해서 변환 정의가 이산 웨이블릿 변환(DWT)이 되는 것은 아님.
DWT는 스케일 및 위치 선택, 필터 구성 등이 다른 별도의 변환임.
선택한 소수의 스케일만으로 원래 신호를 완전히 복원할 수 있다고 보장하지는 않음.

=== 6.5 DNA 인코딩 C 구현
<정의-요약>
이제 비교할 신호와 CWT 적분이 정의되었음.
실제 계산에서는 먼저 아래 함수로 염기 문자를 1.3절의 복소 값으로 바꿈.
C의 `I`는 허수 단위이고 `double complex`는 배정밀도 복소수 형식임.
대소문자를 동일하게 처리하고, A/C/G/T 이외의 문자는 0으로 매핑함.
이 값들로 적분을 계산하는 방법은 다음 장에서 다룸.

```c
static int base_index(char base) {
  switch (base) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

static double complex base_signal(char base) {
  const double complex mapping[4] = {1.0, I, -I, -1.0}; // 이거 나중에 꼭 매핑 바꿔봐야 함.
  int index = base_index(base);
  return index < 0 ? 0.0 : mapping[index];
}
```

== 7. 연속 적분에서 커널 합으로
<연속-적분에서-커널-합으로>
6절의 CWT는 전체 신호와 wavelet을 곱해 적분하는 식임.
이 장의 목표는 그 적분을 #text(weight: "bold")[염기별 기여를 계산하고 더하는 작업]으로 바꾸는 것임.
FFT는 아직 사용하지 않음. 먼저 한 위치에서 계수 하나를 계산하는 방법만 다룸.

=== 7.1 두 염기에서 시작해 적분 나누기
<적분을-염기-구간으로-나눕니다>
예제로 `AC`를 사용함. 2절에서 이 신호를 다음과 같이 정했음.
- A 구간: $- 1 \/ 2 lt.eq t < 1 \/ 2$에서 값이 1임.
- C 구간: $1 \/ 2 lt.eq t < 3 \/ 2$에서 값이 $i$임.
- 나머지 구간: 값이 0임.

스케일 $a$와 비교 위치 $b$를 먼저 하나씩 고정한다고 하자.
계산 순서는 #text(weight: "bold")[A 구간의 기여 계산 → C 구간의 기여 계산 → 두 결과 더하기]임.
구간 안에서 염기 값은 변하지 않으므로 적분 밖으로 꺼낼 수 있음.

$
  W_(upright("AC")) \( a \, b \)
  = 1 / sqrt(a) \[ integral_(- 1 \/ 2)^(1 \/ 2) overline(psi \( \( t - b \) \/ a \)) d t
  + i integral_(1 \/ 2)^(3 \/ 2) overline(psi \( \( t - b \) \/ a \)) d t \] .
$

이 예제에서 적분하는 것은 wavelet이고, 앞에 곱하는 1과 $i$가 염기 값임.
염기가 많아져도 같은 작업을 각 구간에서 반복하면 됨.
염기 $n$의 구간은 $\[ n - 1 \/ 2 \, n + 1 \/ 2 \)$이므로 일반식은

$
  W_x \( a \, b \) = 1 / sqrt(a) sum_(n = 0)^(L - 1) x \[ n \] integral_(n - 1 \/ 2)^(n + 1 \/ 2) overline(psi \( \( t - b \) \/ a \)) thin d t .
$

임. 긴 수식이지만 의미는 #text(weight: "bold")[각 염기 값 × 그 염기 구간의 wavelet 적분을 모두 더함]임.

다음으로 실제 위치 $t$를 Mother wavelet 좌표 $u$로 바꿈.
한 번에 식 전체를 바꾸기보다 다음 세 가지를 따로 확인하면 됨.
+ $u = \( t - b \) \/ a$로 놓음. 즉, 중심 $b$를 빼고 스케일 $a$로 나눔.
+ 구간 양 끝도 같은 방식으로 바꿈. 하한은 $\( n - b - 1 \/ 2 \) \/ a$, 상한은 $\( n - b + 1 \/ 2 \) \/ a$임.
+ $t = b + a u$이므로 $d t = a d u$임. 따라서 앞의 계수는 $a \/ sqrt(a) = sqrt(a)$가 됨.

$
  W_x \( a \, b \) = sqrt(a) sum_(n = 0)^(L - 1) x \[ n \] integral_(\( n - b - 1 \/ 2 \) \/ a)^(\( n - b + 1 \/ 2 \) \/ a) overline(psi \( u \)) thin d u .
$

e.g.) $a = 4$, $b = 0$에서 A의 실제 구간 $\[ - 0.5 \, 0.5 \)$는
Mother wavelet 좌표의 $\[ - 0.125 \, 0.125 \)$가 됨. \
좌표와 구간 폭을 같이 바꾼 것일 뿐, 다른 신호나 다른 정규화 규칙을 선택한 것은 아님.

=== 7.2 반복되는 구간 적분을 커널로 저장
<반복되는-구간-적분을-미리-계산합니다>
앞 식의 적분 구간은 염기 위치 $n$과 출력 위치 $b$ 각각이 아니라 #text(weight: "bold")[차이 $n - b$]에 의존함.
이 차이를 $j$라고 부르자.
- $j = 0$: 비교 중심에 있는 염기.
- $j = - 1$: 중심보다 1 bp 왼쪽의 염기.
- $j = 1$: 중심보다 1 bp 오른쪽의 염기.

같은 스케일에서 중심과의 거리가 같으면 wavelet 적분도 같음.
그러므로 매번 적분하지 않고 거리별로 미리 계산해 저장할 수 있음.
그 값을 $k_a \[ j \]$라고 쓰면

$ k_a \[ j \] = sqrt(a) integral_(\( j - 1 \/ 2 \) \/ a)^(\( j + 1 \/ 2 \) \/ a) overline(psi \( u \)) thin d u $

임. 이 계수 배열이 커널(kernel)임.
커널을 준비한 뒤 한 위치 $b$의 CWT 계수를 구하는 순서는 다음과 같음.
+ 염기 $n$의 값 $x \[ n \]$을 읽음.
+ 중심과의 거리 $j = n - b$에 해당하는 커널 값을 읽음.
+ 두 복소수를 곱함.
+ 모든 염기에서 얻은 곱을 더함.

$ #box(stroke: black, inset: 3pt, [$ W_x \( a \, b \) = sum_n x \[ n \] k_a \[ n - b \] $]) $

e.g.) `AC`에서 중심을 한 칸 옮기면 각 염기의 상대 거리만 바뀜.

#figure(
  table(
    columns: 4,
    table.header([비교 위치], [A의 거리와 기여], [C의 거리와 기여], [합계]),
    [$b = 0$], [$j = 0$: $k_a \[ 0 \]$], [$j = 1$: $i k_a \[ 1 \]$], [$k_a \[ 0 \] + i k_a \[ 1 \]$],
    [$b = 1$], [$j = - 1$: $k_a \[ - 1 \]$], [$j = 0$: $i k_a \[ 0 \]$], [$k_a \[ - 1 \] + i k_a \[ 0 \]$],
  ),
  kind: table,
)

따라서 #text(weight: "bold")[중심을 이동할 때 커널을 다시 적분할 필요는 없음].
같은 스케일의 커널을 모든 위치에서 재사용하면 됨.
여기까지는 정의한 구간별 상수 신호에 대한 정확한 식의 변형임.
실제 커널을 유한 배열로 만들고 수치 적분하는 방법은 8절에서 다룸.

=== 7.3 구간 적분과 점 샘플링 근사의 차이
<점-샘플링-근사와-무엇이-다른가요>
커널을 만들 때 염기 구간 전체를 적분하는 이유를 확인함.
e.g.) `A` 하나가 $n = 0$에 있으면 신호는 점 하나가 아니라 폭 1인 사각 펄스임. \
그 CWT는 $W_x \( a \, b \) = k_a \[ - b \]$이며, 중심에서도 wavelet의 점 하나가 아니라 구간 적분을 사용함.

충분히 큰 스케일에서는 Mother wavelet 좌표의 적분 폭 $1 \/ a$가 작아짐.
그 짧은 구간에서 wavelet이 거의 변하지 않으면 “구간 폭 × 중심의 함숫값”으로 적분을 근사할 수 있음.

$ k_a \[ j \] approx sqrt(a) dot.op 1 / a overline(psi \( j \/ a \)) = 1 / sqrt(a) overline(psi \( j \/ a \)) . $

이 근사가 흔히 사용하는 점 샘플링 커널임.
작은 스케일에서는 한 염기 구간 안에서도 wavelet이 빠르게 진동하므로,
가운데 값 하나로 구간 전체의 기여를 대신하기 어려움.
현재 구현은 이 근사를 사용하지 않고 구간 적분을 수치적으로 계산함.

=== 7.4 커널 생성 C 구현
위 설명과 코드의 대응을 먼저 읽으면 반복문의 의미를 확인하기 쉬움.
- 바깥 반복문: 사용할 스케일을 하나씩 선택함.
- `radius`, `width`: 중심에서 어디까지 저장할지, 배열 원소는 몇 개 필요한지 정함. 계산 근거는 8.1절에서 다룸.
- 안쪽 반복문: 거리 $j$별로 구간 적분을 한 번씩 계산함.
- `sqrt(dilation)`: 7.1절의 변수 치환에서 나온 $sqrt(a)$임.

`dilation`은 $a$, `tap - radius`는 $j$에 대응함.
`morlet_integral`의 본문은 8.3절, `alloc`과 `Wavelets`의 정의는 14절에 있음.

```c
static void wavelets_init(Wavelets *wavelets) {
  memset(wavelets, 0, sizeof(*wavelets));
  for (size_t scale = 0; scale < WAVE_COUNT; scale++) {
    double dilation = wave_scales[scale];
    if (!isfinite(dilation) || dilation <= 0.0 || dilation != floor(dilation))
      fail("CWT scales must be positive natural numbers");
    double radius_value = ceil(8.0 * dilation + 0.5);
    if (radius_value > (INT_MAX - 1) / 2) fail("CWT scale exceeds supported kernel size");
    int radius = (int)radius_value;
    int width = 2 * radius + 1;
    size_t capacity = MAX_WAVE_SIZE > 0 ? (size_t)MAX_WAVE_SIZE : 1;
    if (capacity > (size_t)width) capacity = (size_t)width;
    while (capacity < (size_t)width) capacity *= 2;
    wavelets->kernel[scale] = alloc(capacity, sizeof(double complex));
    wavelets->widths[scale] = width;
    if (width > wavelets->max_width) wavelets->max_width = width;
    for (int tap = 0; tap < width; tap++) {
      double displacement = tap - radius;
      wavelets->kernel[scale][tap] = sqrt(dilation) * morlet_integral(
        (displacement - 0.5) / dilation, (displacement + 0.5) / dilation);
    }
  }
}
```

== 8. 유한 커널과 수치 정확도
<유한-커널과-수치-정확도>
7절에서는 커널 값 하나를 구하려면 wavelet을 한 염기 구간에서 적분해야 한다는 것을 확인했음.
이제 컴퓨터로 그 커널을 만드는 방법을 다룸. 결정할 것은 두 가지임.
+ #text(weight: "bold")[어디까지 저장할 것인가?] 멀리 떨어진 작은 꼬리를 잘라 유한 배열로 만듦.
+ #text(weight: "bold")[각 값을 어떻게 구할 것인가?] 남긴 구간에서 몇 개의 함숫값을 계산해 적분을 근사함.

8.1~8.3절은 이 두 작업을 설명함. 꼬리 오차의 상세 유도는 13.4절에서 따로 다룸.

=== 8.1 무한한 꼬리의 절단 범위
<무한한-꼬리를-어디서-자를까요>
wavelet은 중심에서 멀어질수록 가우스 포락선 때문에 작아짐.
e.g.) $u = 8$에서 포락선 $e^(- u^2 \/ 2)$은 $e^(- 32) approx 1.27 times 10^(- 14)$임. \
정확히 0은 아니지만 매우 작은 값임.
컴퓨터에서 무한 배열을 저장할 수 없으므로 Mother wavelet 좌표에서 $\| u \| > 8$인 부분을 버림.

이 범위를 염기 좌표로 다시 읽으면 $\| t - b \| lt.eq 8 a$임.
e.g.) 스케일 $a = 4$이면 중심에서 양쪽 32 bp까지 남김. \
스케일이 커지면 wavelet 창도 넓어지므로 더 많은 주변 염기를 남겨야 함.

다음으로 “32 bp까지 남김”을 배열 길이로 바꿔야 함.
염기는 점이 아니라 폭 1인 구간이라는 점에 주의함.
+ 거리 $j$의 염기 구간은 $\[ j - 1 \/ 2 \, j + 1 \/ 2 \)$임.
+ 중심점이 범위 밖이어도 구간 일부가 남길 범위와 겹칠 수 있음.
+ 따라서 양쪽으로 반 염기 폭을 더 고려하고, 정수 배열 길이를 위해 올림함.

$ r_a = ceil.l 8 a + 1 \/ 2 ceil.r \, #h(2em) M_a = 2 r_a + 1 . $

$r_a$는 배열의 반경, $M_a$는 왼쪽 $r_a$개 + 중심 1개 + 오른쪽 $r_a$개를 합한 길이임.
$ceil.l v ceil.r$는 $v$ 이상인 가장 작은 정수임.
e.g.) $a = 4$이면 $r_a = ceil.l 32.5 ceil.r = 33$이고 $M_a = 67$임. \
즉, 거리 $- 33$부터 $33$까지 저장함. 바깥쪽 여분 구간은 실제 적분 범위와 겹치지 않아 0이 됨.
자연수 스케일에서는 같은 계산을 $r_a = 8 a + 1$, $M_a = 16 a + 3$으로 쓸 수 있음.

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

=== 8.2 작은 꼬리를 버려도 되는 이유
<꼬리-절단-오차를-유도합니다>
꼬리의 한 점이 작다는 것만으로는 충분하지 않음.
멀리 있는 모든 구간의 기여를 합쳐도 작은지 확인해야 함.
인코딩된 염기 값은 $\| x \( t \) \| lt.eq 1$이므로 버린 부분의 영향을 가우스 꼬리 적분으로 제한할 수 있음.
그 결과, 절단 전후 CWT 계수 차이의 크기는 다음 값 이하임.

$
  #box(stroke: black, inset: 3pt, [$ E_(upright(t a i l)) lt.eq frac(N \( 1 + e^(- 18) \), 4) sqrt(a) thin e^(- 32) $])
$

식의 상세 유도는 13.4절에 있음. 여기서는 크기의 의미만 확인하면 됨.
e.g.) $a = 4$에서 약 $4.76 times 10^(- 15)$, $a = 9$에서 약 $7.14 times 10^(- 15)$ 이하임. \
기본 스케일에서는 꼬리 절단의 영향이 매우 작다는 뜻임.
다만 이는 #text(weight: "bold")[꼬리 절단만의 절대 오차 상한]임.
계수가 거의 0인 경우 상대 오차가 작다는 보장은 아니며,
수치 적분·FFT·complex64 저장 오차까지 포함한 총 오차 상한도 아님.

=== 8.3 적분을 계산하는 실제 방법
<적분을-계산하는-실제-방법>
적분은 구간 전체의 기여를 더하는 계산임.
가장 단순하게는 “구간 폭 × 가운데의 높이”로 근사할 수 있음.
e.g.) 함수가 $\[ - 0.125 \, 0.125 \]$에서 일정하게 2이면 적분은 $0.25 times 2 = 0.5$임. \
하지만 wavelet은 구간 안에서도 변하므로 가운데 한 점만으로는 부족할 수 있음.
현재 구현은 #text(weight: "bold")[여러 위치의 값을 읽고 서로 다른 가중치를 곱해 더하는 방법]을 사용함.
복소 함수에서도 실수부와 허수부에 같은 규칙을 적용하면 됨.

먼저 적분할 구간을 준비함.
+ 한 염기에 대한 하한과 상한을 7절의 식으로 구함.
+ 이 구간에서 $\[ - 8 \, 8 \]$과 겹치는 부분만 남김. 겹치지 않으면 0을 반환함.
+ 남은 길이를 $d$라 하고 $m = ceil.l 2 d ceil.r$개로 나눔.
+ 각 작은 구간의 길이 $d \/ m$은 최대 0.5임.

e.g.) $a = 4$, $j = 0$이면 구간은 $\[ - 0.125 \, 0.125 \]$이고 길이는 0.25임. \
따라서 $m = ceil.l 0.5 ceil.r = 1$로 한 구간만 사용함.
$a = 1$, $j = 0$이면 길이가 1이므로 두 구간으로 나눔.

다음으로 작은 구간 하나에서 16개 위치의 값을 계산함.
미리 정해 둔 위치는 $\[ - 1 \, 1 \]$ 기준이므로 실제 구간 $\[ l \, h \]$로 옮겨야 함.
- 중심 $v = \( l + h \) \/ 2$: 어디로 이동할지 정함.
- 반길이 $q = \( h - l \) \/ 2$: 얼마나 줄이거나 늘릴지 정함.
- 기준 위치 $z_j$를 $u_j = v + q z_j$로 옮김.
- 그 위치의 함숫값에 가중치 $w_j$를 곱해 더하고, 구간 크기를 반영하는 $q$를 곱함.

$
  integral_l^h f \( u \) thin d u = q integral_(- 1)^1 f \( v + q z \) thin d z approx q sum_(j = 1)^16 w_j f \( v + q z_j \) .
$

e.g.) 위 스케일 4의 중심 구간에서는 $v = 0$, $q = 0.125$임. \
기준 위치 $z_j = 0.5$를 예로 들면 실제 위치는 $u_j = 0.0625$가 됨.
이는 좌표 이동을 보여주는 예일 뿐이며 실제 16개 위치는 아래 `nodes`의 양수·음수 쌍임.

이 방법을 16점 Gauss-Legendre 구적법(quadrature)이라고 부름.
31차 이하의 다항식에는 정확한 값을 주지만, 일반적인 함수에서 오차가 없는 공식은 아님.
현재처럼 매끄러운 함수를 짧은 구간에서 적분하면 좋은 근사를 얻을 수 있음.
위치와 가중치 자체를 유도하는 방법은 13절에서 다룸.

실제 구현은 대칭인 양의 위치 8개만 저장하고 양쪽 위치를 모두 계산함.
`segments`는 $m$, `middle`은 $v$, `step / 2.0`은 $q$에 대응함.
`cexp(-I * 6.0 * ...)`의 음수 부호는 #text(weight: "bold")[켤레 wavelet을 적분]하기 때문에 들어감.
별도로 다시 켤레를 적용하면 안 됨.

```c
static double complex morlet_integral(double lower, double upper) {
  static const double nodes[8] = {
    0.09501250983763744, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
    0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499
  };
  static const double weights[8] = {
    0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541
  };
  lower = fmax(lower, -8.0);
  upper = fmin(upper, 8.0);
  if (lower >= upper) return 0.0;
  double normalization = 1.0 / sqrt(sqrt(acos(-1.0)) *
                                   (1.0 + exp(-36.0) - 2.0 * exp(-27.0)));
  int segments = (int)ceil((upper - lower) * 2.0);
  double step = (upper - lower) / segments;
  double complex integral = 0.0;
  for (int segment = 0; segment < segments; segment++) {
    double middle = lower + (segment + 0.5) * step;
    for (int node = 0; node < 8; node++) {
      double left = middle - nodes[node] * step / 2.0;
      double right = middle + nodes[node] * step / 2.0;
      integral += weights[node] * step / 2.0 * (
        exp(-0.5 * left * left) * (cexp(-I * 6.0 * left) - exp(-18.0)) +
        exp(-0.5 * right * right) * (cexp(-I * 6.0 * right) - exp(-18.0)));
    }
  }
  return normalization * integral;
}
```

=== 8.4 계산 정밀도와 출력 형식
<구현-요약>
7.4절의 `wavelets_init`이 거리별 구간을 정하고, 8.3절의 `morlet_integral`이 그 구간을 적분함.
이 결과에 $sqrt(a)$를 곱하면 커널 값이 완성됨.
컴퓨터 계산에서 발생하는 근사는 서로 구분해야 함.
- 꼬리 절단: 먼 구간을 버리는 근사.
- 수치 적분: 구간 전체를 유한한 평가 위치로 대신하는 근사.
- 부동소수점 연산: 계산 중 유한한 정밀도로 반올림하는 오차.
- 출력 형식: 배정밀도 계산값을 complex64로 바꾸는 반올림.

complex64는 실수부와 허수부를 각각 32비트 부동소수점으로 저장하는 형식임.
따라서 수식의 정의는 명확하더라도 파일에 기록된 숫자가 오차 없는 정확 산술의 결과는 아님.

== 9. 커널 합을 FFT로 계산하기
<커널-합을-fft로-계산하기>
7~8절까지의 결과는 “커널을 준비하고 위치마다 곱한 값들을 더함”임.
이대로 구현해도 CWT를 계산할 수 있음. 문제는 긴 서열의 모든 위치와 스케일에서 반복하면 계산량이 커진다는 점임.
9절에서는 #text(weight: "bold")[같은 결과를 더 빠르게 얻는 계산 순서]를 다룸.
새로운 CWT 정의를 만드는 것이 아님.

먼저 필요한 작업을 요약하면 #text(weight: "bold")[커널 뒤집기 → 0 채우기 → FFT → 주파수별 곱 → 역 FFT → 필요한 위치 읽기]임.
아래에서는 각 작업이 왜 필요한지 하나씩 확인함.

=== 9.1 합성곱에 맞게 커널 방향 변경
<커널-방향을-먼저-맞춥니다>
합성곱(convolution)은 두 배열에서 값의 곱을 만들어 특정 규칙으로 더하는 계산임.
작은 숫자 배열로 규칙부터 확인함. 여기서 숫자는 연습용이며 실제 Morlet 커널은 아님.
입력을 $X = \[ 1, 2, 3 \]$, 커널을 $H = \[ 10, 20 \]$이라 하면 배열 밖의 값은 0임.

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

즉, 합성곱 결과는 $Y = \[ 10, 40, 70, 60 \]$임.
각 출력에서 사용하는 입력 인덱스와 커널 인덱스의 합이 $m$이 된다는 점에 주목함.
이를 일반적인 신호 위치 $b$에 대해 쓰면

$ \( x * h \) \[ b \] = sum_n x \[ n \] h \[ b - n \] $

임. 입력 인덱스가 $n$이므로 커널 인덱스는 $b - n$임.
그러나 CWT의 커널 인덱스는 $n - b$로 부호가 반대임.
이 차이만 해결하면 이미 알려진 합성곱 계산법을 사용할 수 있음.

$ h_a \[ j \] = k_a \[ - j \] $

로 커널의 좌우를 뒤집으면 $h_a \[ b - n \] = k_a \[ n - b \]$가 됨.
따라서 CWT와 합성곱이 같은 계산이 됨.

#text(weight: "bold")[켤레와 배열 뒤집기는 서로 다른 작업]임.
- 켤레: 허수부의 부호를 바꾸는 작업. CWT 정의에서 이미 적용함.
- 배열 뒤집기: 위치 $j$와 $- j$의 계수를 교환하는 작업. 합성곱 형식에 맞추기 위해 적용함.

두 번 켤레를 적용하거나 커널 방향을 그대로 두면 다른 변환을 계산하게 됨.

=== 9.2 음수 인덱스를 배열에 저장하는 방법
<음수-인덱스를-배열에-저장하는-방법>
수학에서는 중심을 0으로 두어 커널의 왼쪽을 음수 인덱스로 쓸 수 있음.
그러나 C 배열의 인덱스는 0부터 시작함.
따라서 #text(weight: "bold")[수학적 거리와 배열 인덱스를 구분]해야 함.
e.g.) 반경 $r = 1$인 커널은 다음과 같이 저장됨.

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

반경이 $r$일 때도 같음.
거리 $- r, dots.h, r$를 배열 인덱스 $0, dots.h, 2 r$에 대응시키므로
$K \[ q \] = k_a \[ q - r \]$임. 이 배열을 뒤집으면

$ H \[ q \] = K \[ 2 r - q \] = k_a \[ r - q \] $

임. 중심이 배열의 $r$번에 있으므로 합성곱 결과를 읽을 때도 그만큼의 이동을 반영해야 함.
입력 배열이 contig의 위치 $s$부터 시작한다면, 원하는 위치 $b$의 계수는
#text(weight: "bold")[입력 안의 위치 $b - s$ + 커널 중심의 배열 위치 $r$]에서 읽음.

$ #box(stroke: black, inset: 3pt, [$ m = b - s + r $]) $

e.g.) 입력을 위치 100부터 읽었고, 위치 105의 계수가 필요하며, 커널 반경이 9이면
$m = 105 - 100 + 9 = 14$임. \
즉, 합성곱 배열의 인덱스 14를 읽어야 함. 5를 그대로 읽으면 다른 위치의 결과가 됨.

식으로 확인하려면 $X \[ p \] = x \[ s + p \]$를 대입함.
$p$는 입력 배열 안의 상대 위치이며, 일반 합성곱은

$ Y \[ m \] = sum_p X \[ p \] H \[ m - p \] = sum_p x \[ s + p \] k_a \[ r - m + p \] . $

임. 원하는 커널 인덱스 $s + p - b$와 같아지도록 $r - m + p = s + p - b$를 풀면 위의 $m$을 얻음.
코드의 `width - 1 - width/2`는 홀수 길이 $2 r + 1$에서 정확히 $r$임.

=== 9.3 DFT로 신호를 회전 패턴에 분해
<dft는-신호를-회전-패턴들로-분해합니다>
FFT를 이해하기 전에 DFT(discrete Fourier transform; 이산 푸리에 변환)를 정의함.
DFT는 위치별 값을 #text(weight: "bold")[어떤 주파수의 진동이 얼마나 들어 있는가]라는 표현으로 바꾸는 계산임.
입력도 길이 $P$, 출력도 길이 $P$인 복소 배열이며, 역변환으로 원래 배열을 되찾을 수 있음.
이 주파수 좌표로 바꾸면 합성곱을 더 간단히 계산할 수 있으므로 사용함.

e.g.) $X = \[ 1, 1, 1, 1 \]$은 위치에 따라 변하지 않는 신호임. \
DFT 결과는 $\[ 4, 0, 0, 0 \]$이며 상수 성분에만 값이 생김.
다른 진동과 비교하면 양쪽 기여가 상쇄되기 때문임.

길이 $P$인 배열의 정방향 변환식은

$ hat(X) \[ k \] = sum_(n = 0)^(P - 1) X \[ n \] e^(- 2 pi i k n \/ P) \, #h(2em) k = 0 \, dots.h \, P - 1 . $

각 $k$에 대해 다른 회전 속도로 배열을 비교한 값임.
$k$는 길이 $P$ 안에서의 회전 횟수에 해당하는 주파수 인덱스임.
간격이 1 bp이면 $k \/ P$ cycles/bp에 대응하되 높은 인덱스는 음의 주파수로도 해석됨.
e.g.) $k = P - 1$은 $- 1 \/ P$ cycles/bp와 같은 샘플 회전을 나타냄.
역변환은

$ X \[ n \] = 1 / P sum_(k = 0)^(P - 1) hat(X) \[ k \] e^(2 pi i k n \/ P) $

임. 앞의 $1 \/ P$는 크기를 원래대로 되돌리는 계수임.
상수 신호 예제에서도 값 4를 길이 4로 나누면 원래 값 1을 얻음.
역변환이 성립하는 상세 유도는 11.3절에서 다룸.

=== 9.4 주파수별 곱셈이 합성곱이 되는 이유
<주파수별-곱셈이-합성곱이-되는-이유>
사용할 성질은 다음과 같음.
+ 두 배열을 같은 길이로 맞추고 각각 DFT함.
+ 같은 주파수 인덱스끼리 복소수를 곱함.
+ 그 결과를 역변환하면 두 배열의 합성곱을 얻음.

다만 DFT는 배열의 끝과 처음이 이어지는 주기 신호로 취급함.
따라서 별도 조치가 없으면 #text(weight: "bold")[원형 합성곱(circular convolution)]이 계산됨.
이는 끝을 넘어간 결과가 처음으로 돌아오는 계산임.
성질의 수식 유도는 11.4절에서 다룸.

e.g.) 9.1절에서 얻은 결과 $\[ 10, 40, 70, 60 \]$을 길이 3에 담으려 하면,
마지막 60이 처음으로 돌아가 $\[ 70, 40, 70 \]$이 됨. \
우리가 원하는 CWT에는 이처럼 서열 끝의 신호가 처음으로 돌아오는 효과가 없어야 함.

우리가 원하는 것은 끝과 처음이 연결되지 않는 #text(weight: "bold")[선형 합성곱(linear convolution)]임.
이를 위해 입력과 커널 뒤에 0을 채워 계산 공간을 넓힘.
입력 길이가 $S$, 커널 길이가 $M$이면 가능한 마지막 출력 인덱스는
$\( S - 1 \) + \( M - 1 \)$임. 0번부터 세므로 필요한 출력 길이는 $S + M - 1$임.
따라서 계산 길이 $P$를

$ P gt.eq S + M - 1 $

이 되도록 두 배열 뒤에 0을 채우면, 0이 아닌 결과가 반대편으로 돌아 겹치지 않음.
이를 0 채우기(zero-padding)라고 부름.
현재 구현은 최대 커널 길이를 기준으로 이 조건을 만족하는 가장 작은 2의 거듭제곱 $P$를 선택함.

e.g.) 앞 예제는 $S = 3$, $M = 2$이므로 길이 4가 필요함. \
$P = 4$로 잡고 $X = \[ 1, 2, 3, 0 \]$, $H = \[ 10, 20, 0, 0 \]$으로 준비함.
이 두 배열을 변환·곱셈·역변환하면 원하는 $\[ 10, 40, 70, 60 \]$을 얻음.
$S = 4$, $M = 3$인 다른 입력에서는 길이 6이 필요하므로 $P = 8$을 선택함.

`rfft.h`의 역변환은 위 식의 $1 \/ P$를 적용하지 않으므로 결과를 #text(weight: "bold")[FFT 길이 $P$로] 나눔.
$S + M - 1$과 $P$는 일반적으로 다름. 둘을 혼동하면 진폭이 틀려짐.

FFT 길이 선택 구현:

```c
static size_t convolution_fft_size(size_t signal_length, size_t kernel_length) {
  if (signal_length > SIZE_MAX - (kernel_length - 1)) fail("CWT FFT size overflow");
  size_t required = signal_length + kernel_length - 1;
  size_t size = 1;
  while (size < required) {
    if (size > SIZE_MAX / 2) fail("CWT FFT size overflow");
    size <<= 1;
  }
  return size;
}
```

=== 9.5 FFT는 DFT를 빠르게 계산하는 알고리즘
<fft는-새로운-수학적-변환이-아닙니다>
DFT 식을 그대로 계산하면 약 $P^2$개의 곱셈이 필요함.
FFT(fast Fourier transform; 고속 푸리에 변환)는 같은 합을 짝수 위치와 홀수 위치로 나눠 재사용하는 계산 방법임.
길이 $P$를 절반씩 나눠 계산하므로 현재 구현에서는 2의 거듭제곱 길이를 사용함.
이 분할을 반복하면 계산량이 대략 $P log_2 P$로 줄어듦.
e.g.) $P = 1024$이면 $P^2$는 약 백만인 반면 $P log_2 P$는 약 만임. \
실제 연산 횟수는 구현에 따라 달라지지만 긴 배열에서 속도 차이가 커지는 이유를 보여줌.
변환의 정의는 그대로이며 계산 순서만 달라짐. 분할식은 11.5절에 있음.

CWT 계산의 전체 흐름을 다시 연결하면 다음과 같음.
+ 7~8절에서 만든 커널을 좌우로 뒤집음. 이때 다시 켤레를 적용하지 않음.
+ 입력과 커널의 합성곱 전체를 담을 길이 $P$를 선택하고 두 배열 뒤에 0을 채움.
+ 입력과 커널을 각각 FFT해 위치 좌표에서 주파수 좌표로 옮김.
+ 같은 주파수 인덱스의 복소수끼리 곱함.
+ 역 FFT로 위치 좌표로 돌아옴. 현재 라이브러리 결과를 $P$로 나눠 크기를 맞춤.
+ 원하는 위치 $b$마다 합성곱 배열의 $b - s + r$번 값을 읽음.

같은 입력 청크의 FFT는 여러 스케일에서 재사용할 수 있음.
스케일마다 바뀌는 것은 커널과 그 FFT임.

=== 9.6 청크 경계와 contig 경계의 구분
<청크-경계와-contig-경계는-다릅니다>
메모리 사용량을 제한하기 위해 프로그램은 출력 위치를 1024개씩 나눠 처리함.
이 묶음을 청크(chunk)라고 부름.
청크는 계산을 나누는 단위일 뿐 새로운 contig가 아님.
청크의 맨 앞이나 뒤 위치를 계산할 때도 이웃 청크의 염기가 필요함.
따라서 출력할 구간만 읽지 않고 양쪽으로 최대 커널 반경만큼 더 읽음.

e.g.) 출력할 위치가 $1000, dots.h, 2023$이고 최대 반경이 73이면,
입력은 위치 927부터 2096까지 필요함. \
추가로 읽은 부분은 주변 문맥이며, 이번 청크에서는 원래 요청한 1024개 위치의 결과만 저장함.
contig 밖으로 나간 부분에는 실제 염기가 없으므로 0을 사용함.

청크 시작을 $s_0$, 출력 개수를 $B$, 최대 반경을 $R$이라 하면 입력 문맥은

$ \[ s_0 - R \, #h(0em) s_0 + B + R \) $

과 contig 구간의 교집합임.
이 문맥이 있어야 청크 경계에서도 전체 신호를 한 번에 계산한 것과 같은 결과를 얻음.
- 청크 경계: 계산 편의를 위해 나눈 경계. 그 밖의 신호도 존재하므로 주변 염기를 읽음.
- contig 경계: 입력 서열 자체의 경계. 신호 정의에 따라 그 밖은 0으로 처리함.

실제 CWT 계산 구현은 다음과 같음.
`context_start`는 $s$, `position`은 $b$, `fft_size`는 $P$에 대응함.
`alloc`이 0으로 초기화하므로 입력을 채우고 남은 부분은 자동으로 zero-padding이 됨.
`features`에는 위치별로 스케일 순서의 실수부·허수부를 번갈아 저장함.

```c
static void cwt_extract(const Wavelets *wavelets, const char *sequence, int length,
                        int start, int count, double *features) {
  if (count <= 0) return;
  memset(features, 0, (size_t)count * CWT_CHANNELS * sizeof(double));
  if (start >= length || (long)start + count <= 0) return;

  long context_start = (long)start - wavelets->max_width / 2;
  long context_end = (long)start + count + wavelets->max_width / 2;
  if (context_start < 0) context_start = 0;
  if (context_end > length) context_end = length;
  size_t signal_length = (size_t)(context_end - context_start);
  size_t fft_size = convolution_fft_size(signal_length, (size_t)wavelets->max_width);
  double complex *signal_fft = alloc(fft_size, sizeof(*signal_fft));
  double complex *kernel_fft = alloc(fft_size, sizeof(*kernel_fft));

  for (size_t index = 0; index < signal_length; index++)
    signal_fft[index] = base_signal(sequence[context_start + (long)index]);
  fft_transform(signal_fft, fft_size, false);

  for (size_t scale = 0; scale < WAVE_COUNT; scale++) {
    int width = wavelets->widths[scale];
    memset(kernel_fft, 0, fft_size * sizeof(*kernel_fft));
    for (int tap = 0; tap < width; tap++)
      kernel_fft[width - 1 - tap] = wavelets->kernel[scale][tap];
    fft_transform(kernel_fft, fft_size, false);
    for (size_t index = 0; index < fft_size; index++) kernel_fft[index] *= signal_fft[index];
    fft_transform(kernel_fft, fft_size, true);

    for (int offset = 0; offset < count; offset++) {
      long position = (long)start + offset;
      if (position < 0 || position >= length) continue;
      size_t convolution_index = (size_t)(position - context_start) + width - 1 - width / 2;
      double complex coefficient = kernel_fft[convolution_index] / (double)fft_size;
      features[(size_t)offset * CWT_CHANNELS + 2 * scale] = creal(coefficient);
      features[(size_t)offset * CWT_CHANNELS + 2 * scale + 1] = cimag(coefficient);
    }
  }

  free(signal_fft);
  free(kernel_fft);
}
```

== 10. 출력 행렬과 결과 해석
<설정과-출력>
앞 장에서는 위치별 CWT 계수를 계산했음.
이제 그 계수가 파일의 어디에 저장되는지 확인하고, 행렬을 그림으로 읽는 방법을 다룸.

=== 10.1 스케일 설정과 출력 형식
`wave_scales`는 `CWT_SCALES`로 설정하며 기본값은 `4,5,6,7,8,9`임.
스케일은 양의 자연수여야 하며 개수·순서는 자유롭게 지정할 수 있고 중복도 허용함.
커널 크기는 8.1절의 반경으로 정해지며, 지원 범위와 사용 가능한 메모리 안에서 스케일을 선택해야 함.

e.g.) 다음과 같이 설정할 수 있음. 출력 디렉터리는 아직 존재하지 않는 경로여야 함.

```sh
cc -O3 -std=c11 -Wall -Wextra '-DCWT_SCALES=1,2,8,32' \
  anno_cwt.c -lz -lm -o anno_cwt_custom
./anno_cwt_custom genome.fa new_output_directory
```

출력 파일:
- `contigs.tsv`: 스케일 순서, contig 이름과 길이, 각 가닥의 행 시작 위치.
- `matrix.bin`: little-endian complex64 형식의 복소 CWT 계수.

=== 10.2 행과 열의 의미
스케일 개수를 $J$, 모든 contig 길이의 합을 $L_(upright("total"))$이라 하면
저장 행렬의 크기는 $2 L_(upright("total")) times J$임.
- 행(row): contig와 가닥을 포함한 염기 위치.
- 열(column): 설정한 순서의 스케일.
- 원소: 해당 위치와 스케일의 복소 CWT 계수.

contig 하나에 대해 양의 가닥의 $L$개 행 다음에 역상보 서열의 $L$개 행을 저장함.
`plus_offset`과 `minus_offset`은 각 가닥 블록의 #text(weight: "bold")[시작 행 인덱스]이며 바이트 위치가 아님.
가닥 블록의 시작 행을 $o$, 스케일 순서를 $a_0, dots.h, a_(J - 1)$이라 하면

$ M \[ o + b \, j \] = W_x \( a_j \, b \) . $

e.g.) 기본 설정의 첫 번째 열은 스케일 4, 마지막 열은 스케일 9임. \
열 번호 자체가 스케일 값은 아니므로 반드시 `scales` 필드를 함께 읽어야 함.
전체 파일 크기는 $2 L_(upright("total")) dot.op J dot.op 8$ 바이트임.

=== 10.3 DNA CWT 행렬과 viridis 그림 읽기
앞 절에서 원소 하나가 특정 위치와 스케일의 복소 CWT 계수임을 확인했음.
이제 #text(weight: "bold")[계수를 파워로 바꾸기 → 위치와 스케일에 배치하기 → 밝은 패턴 읽기] 순서로 그림을 해석함.

==== 10.3.1 파워를 색으로 표시
원소 하나가 $W = R + i I$이면 파워는 실수부·허수부의 제곱합임.

$ P \[ b \, j \] = \| W_x \( a_j \, b \) \|^2 = R^2 + I^2 . $

e.g.) $W = 1 + i$이면 $P = 2$, $W = - 1 - i$여도 $P = 2$임. \
두 계수의 위상은 다르지만 파워 색은 같음.
현재 viridis 패널은 큰 값의 차이를 압축하기 위해 $D = ln \( 1 + P \)$를 색으로 표시함.
여러 위치를 한 픽셀로 줄이는 경우에는 파워를 먼저 평균한 뒤 로그를 적용함.
색은 #text(weight: "bold")[어두운 보라 → 청록 → 노랑] 순서로 값이 커짐.
같은 색이라도 그림마다 표시 범위가 다를 수 있으므로 그림 간 비교에는 컬러바의 숫자를 사용함.
상한에서 포화된 색만으로는 그보다 큰 값들의 차이를 알 수 없음.

==== 10.3.2 위치와 스케일을 축으로 읽기
위치별·스케일별 파워를 표시한 그림을 스칼로그램(scalogram)이라고 부름.
저장 행렬에서는 행이 위치, 열이 스케일이지만, 그림에서는 보통 전치하여 위치를 가로축에 놓음.
이때 가로축을 따라가면 같은 스케일의 응답이 서열을 따라 어떻게 달라지는지 볼 수 있음.
세로축을 따라가면 같은 위치가 어떤 스케일에 강하게 반응하는지 볼 수 있음.

스케일은 주기 자체가 아님. 현재 Morlet의 대표 주기는

$ T approx 2 pi a \/ 6 upright(" bp") $

임. 유도는 12.2절에서 다룸.
e.g.) 스케일 4의 대표 주기는 약 4.19 bp, 스케일 9는 약 9.42 bp임. \
기본 스케일 4~9는 이 범위의 주기 성분을 비교하는 설정임.
각 wavelet은 주변 주파수에도 반응하므로 밝은 스케일 하나가 정확한 반복 길이 하나를 결정하지는 않음.

c.f.) 여러 구간을 정렬한 그림에서는 세로축이 스케일이 아니라 구간 순위일 수도 있음. \
아래 설명은 위치와 스케일을 두 축으로 놓은 스칼로그램을 기준으로 함.

==== 10.3.3 밝은 띠와 피크의 의미
- #text(weight: "bold")[한 스케일에서 밝은 띠가 이어짐]: 그 스케일이 잘 반응하는 진동 성분이 일정 구간에서 지속될 가능성이 있음.
- #text(weight: "bold")[한 위치에서 여러 스케일이 동시에 밝음]: 급격한 변화나 여러 주파수 성분의 혼합일 수 있음.
- #text(weight: "bold")[위치를 따라 밝은 띠의 스케일이 달라짐]: 구간별로 지배적인 진동 성분이 바뀔 가능성이 있음.
- #text(weight: "bold")[어두운 구간]: 선택한 인코딩과 스케일에서 반응이 약하다는 뜻임.

계수 하나는 중심 염기만이 아니라 주변 구간을 함께 반영함.
따라서 1 bp마다 출력한다고 변화의 위치를 1 bp 정확도로 구분하는 것은 아님.
큰 스케일에서는 분석 창과 피크가 더 넓어질 수 있음.
또한 에너지 정규화는 비교 wavelet에 적용한 것이므로,
동일한 진폭의 입력 패턴이 모든 스케일에서 같은 파워를 갖는다고 보장하지는 않음.

큰 파워는 해당 wavelet에 대한 강한 응답이지 유전자나 생물학적 활성의 직접적인 지표는 아님.
반복 서열, 염기 조성 변화, 모호한 염기의 0 처리, 입력 경계에도 반응할 수 있음.
패턴을 판단할 때는 원래 서열과 함께 확인함. 가닥 좌표와 입력 경계의 영향은 다음 절에서 다룸.

=== 10.4 가닥 좌표와 행렬 저장 구현
출력 위치 $b$는 각 가닥의 서열에서 0부터 시작함.
원래 contig 길이를 $L$, 1부터 시작하는 유전체 좌표를 $g$라 하면
- 양의 가닥: $g = b + 1$.
- 음의 가닥: $g = L - b$.

음의 가닥은 출력 행만 뒤집은 것이 아니라 #text(weight: "bold")[역상보 서열을 별도로 인코딩하고 변환]한 결과임.
유전체 좌표의 동일 구간을 비교할 때는 가닥별 행 시작 위치와 방향을 모두 확인해야 함.
contig 끝에서는 wavelet 일부가 0 확장 영역과 겹치므로 내부 위치와 다른 응답이 나올 수 있음.
큰 스케일일수록 이 경계 효과가 나타나는 범위도 넓어짐.
전체 contig의 계수를 계산한 뒤 관심 구간만 잘라 표시하는 것은 새로운 입력 경계를 만들지 않음.
반면 FASTA 자체를 짧게 잘라 다시 계산하면 잘린 끝이 입력 경계가 되므로 결과가 달라질 수 있음.

저장 구현은 다음과 같음.
배정밀도 `features`를 32비트 실수 배열 `values`로 바꾼 뒤,
실수부·허수부 쌍을 complex64의 한 원소로 저장함.
`read_fasta`, `create_in`, `write_floats`는 14절에 포함함.

```c
static void export_cwt(const Wavelets *wavelets, const char *fasta, const char *directory) {
  int count;
  Contig *contigs = read_fasta(fasta, &count);
  if (mkdir(directory, 0777)) fail("Cannot create directory %s", directory);
  FILE *matrix = create_in(directory, "matrix.bin");
  FILE *index = create_in(directory, "contigs.tsv");
  fputs("# anno-cwt-v2\n# scales\t", index);
  for (size_t scale = 0; scale < WAVE_COUNT; scale++) fprintf(index, "%s%.17g", scale ? "," : "", wave_scales[scale]);
  fputs("\nseqid\tlength\tplus_offset\tminus_offset\n", index);
  double *features = alloc(1024 * CWT_CHANNELS, sizeof(*features));
  float *values = alloc(1024 * CWT_CHANNELS, sizeof(*values));
  size_t rows = 0;
  for (int contig_index = 0; contig_index < count; contig_index++) {
    Contig *contig = contigs + contig_index;
    fprintf(index, "%s\t%d\t%zu\t%zu\n", contig->name, contig->length, rows, rows + contig->length);
    char *reverse = alloc((size_t)contig->length + 1, 1);
    for (int position = 0; position < contig->length; position++) {
      int base = base_index(contig->seq[contig->length - 1 - position]);
      reverse[position] = base < 0 ? 'N' : "TGCA"[base];
    }
    for (int strand = 0; strand < 2; strand++) {
      const char *bases = strand ? reverse : contig->seq;
      for (int start = 0; start < contig->length; start += 1024) {
        int length = contig->length - start < 1024 ? contig->length - start : 1024;
        cwt_extract(wavelets, bases, contig->length, start, length, features);
        for (size_t value = 0; value < (size_t)length * CWT_CHANNELS; value++) values[value] = (float)features[value];
        write_floats(matrix, values, (size_t)length * CWT_CHANNELS);
      }
    }
    free(reverse);
    rows += 2 * (size_t)contig->length;
  }
  if (fclose(matrix)) fail("Cannot finish CWT matrix");
  if (fclose(index)) fail("Cannot finish CWT index");
  for (int contig_index = 0; contig_index < count; contig_index++) {
    free(contigs[contig_index].name);
    free(contigs[contig_index].seq);
  }
  free(contigs);
  free(features);
  free(values);
  fprintf(stderr, "Saved %zu rows x %zu complex64 scales\n", rows, WAVE_COUNT);
}
```

== 11. 보충: 가우스 적분과 푸리에 변환
<가우스-적분을-처음부터-유도하기>
이 장에서는 4~5절의 가우스 적분과 9절의 DFT·FFT 성질을 유도함.
본문에서는 계산 방법을 먼저 설명했으며, 여기서는 그 방법이 성립하는 근거를 확인함.
가우스 적분을 상수와 함수의 식으로 표현하려면 평면의 극좌표라는 추가 개념이 필요함.
이를 사용하지 않고 수치 적분으로 정규화 계수를 근사할 수도 있지만,
현재 구현의 $N$과 $e^(- 18)$이라는 식은 아래 계산에서 나옴.

=== 11.1 가우스 곡선 아래 면적
<가우스-곡선-아래-면적>
$J = integral_(bb(R)) e^(- u^2 \/ 2) thin d u$라고 놓음.
같은 적분을 두 번 곱하고 두 적분의 변수를 $x$, $y$로 구분하면

$ J^2 = integral_(bb(R)) integral_(bb(R)) e^(- \( x^2 + y^2 \) \/ 2) thin d x thin d y . $

평면 전체에서 높이 $e^(- \( x^2 + y^2 \) \/ 2)$를 면적에 곱해 더하는 적분임.
높이가 원점으로부터의 거리에만 의존하므로 극좌표로 바꾸면 계산이 쉬워짐.
- $x = r cos theta$, $y = r sin theta$로 놓으면 $x^2 + y^2 = r^2$임.
- $r$은 원점으로부터의 거리, $theta$는 각도임.
- 반지름 폭 $d r$, 각도 폭 $d theta$인 작은 부채꼴의 면적은
  호 길이 $r thin d theta$와 폭 $d r$의 곱인 $r thin d r thin d theta$임.

따라서

$
  J^2 = integral_0^(2 pi) integral_0^oo e^(- r^2 \/ 2) r thin d r thin d theta = 2 pi \[ - e^(- r^2 \/ 2) \]_0^oo = 2 pi .
$

$J > 0$이므로 $J = sqrt(2 pi)$임.
이 문서에서는 다중적분과 좌표변환의 일반 이론 대신 작은 면적을 더하는 방식으로 이해함.

=== 11.2 가우스 함수와 진동을 곱한 적분
<가우스-함수와-진동을-곱한-적분>
각주파수 $s$에 따라 적분값이 어떻게 달라지는지 보기 위해 다음 함수를 정의함.

$ F \( s \) = integral_(bb(R)) e^(- u^2 \/ 2) e^(i s u) thin d u . $

가우스 함수는 함수 자체도, $\| u \|$를 곱한 함수도 적분 가능할 정도로 빠르게 감소함.
이 성질을 이용하면 $s$에 대한 미분을 적분 안으로 옮길 수 있음.
$e^(i s u)$를 $s$로 미분하면 $i u e^(i s u)$이므로

$ F' \( s \) = i integral_(bb(R)) u e^(- u^2 \/ 2) e^(i s u) thin d u . $

$g \( u \) = e^(- u^2 \/ 2)$이면 $g' \( u \) = - u g \( u \)$이므로 부분적분을 적용함.

$
  F' \( s \) & = - i integral_(bb(R)) g' \( u \) e^(i s u) thin d u \
             & = - i \[ g \( u \) e^(i s u) \]_(- oo)^oo + i integral_(bb(R)) g \( u \) \( i s \) e^(i s u) thin d u \
             & = - s F \( s \) .
$

경계항은 $\| e^(i s u) \| = 1$이고 $g \( u \) arrow.r 0$이므로 0임.
결국 $F' \( s \) = - s F \( s \)$이라는 관계를 얻음.
곱의 미분법을 사용하면

$ frac(d, d s) \[ e^(s^2 \/ 2) F \( s \) \] = e^(s^2 \/ 2) \[ s F \( s \) + F' \( s \) \] = 0 . $

따라서 $e^(s^2 \/ 2) F \( s \)$의 값은 상수임.
$F \( 0 \) = sqrt(2 pi)$를 대입하면

$ #box(stroke: black, inset: 3pt, [$ F \( s \) = sqrt(2 pi) e^(- s^2 \/ 2) $]) . $

실수부를 취하면 코사인 적분 항등식을 얻음.
허수부를 취하면 사인 적분이 0임을 알 수 있음.
가우스 함수는 짝함수이고 사인은 홀함수이므로 양쪽 구간의 기여가 상쇄된다는 설명과도 일치함.
또한 $v = sqrt(2) u$로 치환하면

$ integral_(bb(R)) e^(- u^2) e^(i s u) thin d u = sqrt(pi) e^(- s^2 \/ 4) $

임. 여기에 $s = 6$을 대입한 것이 5절의 $sqrt(pi) e^(- 9)$임.

=== 11.3 DFT 역변환에서 원래 값만 남는 이유
9.3절의 역변환을 확인하는 핵심은 회전하는 복소수들의 등비수열 합임.

$
  sum_(k = 0)^(P - 1) e^(2 pi i k \( n - m \) \/ P) = cases(delim: "{", P \, & n = m med \( mod med P \) \,, 0 \, & upright("그 밖의 경우") .)
$

$n = m med \( mod med P \)$는 두 인덱스의 차이가 $P$의 정수배라는 뜻임.
회전비가 1이면 항 $P$개의 합은 $P$임.
그 외에는 등비수열 공식의 분자가 $1 - e^(2 pi i \( n - m \)) = 0$이므로 합이 0임.
9.3절의 정방향 식을 역변환에 대입하면 이 성질 때문에 원래 위치의 값만 남음.

=== 11.4 주파수별 곱이 원형 합성곱인 이유
두 DFT를 곱해 역변환에 대입하고 합의 순서를 바꾸면

$
  1 / P sum_k hat(X) \[ k \] hat(H) \[ k \] e^(2 pi i k m \/ P) = sum_p sum_q X \[ p \] H \[ q \] \( 1 / P sum_k e^(2 pi i k \( m - p - q \) \/ P) \) .
$

11.3절의 성질에 의해 안쪽 합을 $P$로 나눈 값은 $m = p + q med \( mod med P \)$일 때만 1이고 나머지는 0임.
그러므로 입력 인덱스와 커널 인덱스의 합이 출력 인덱스와 같은 항만 남음.
길이 $P$마다 인덱스가 돌아온다는 점을 제외하면 9.1절의 합성곱 규칙과 같음.
9.4절의 0 채우기는 이 돌아오는 항들이 원래 결과에 겹치지 않게 하는 작업임.

=== 11.5 FFT에서 길이를 절반으로 나누는 이유
$zeta = e^(- 2 pi i \/ P)$로 놓고 입력의 짝수 위치와 홀수 위치를 분리하면

$ hat(X) \[ k \] = sum_j X \[ 2 j \] \( zeta^2 \)^(k j) + zeta^k sum_j X \[ 2 j + 1 \] \( zeta^2 \)^(k j) . $

오른쪽의 두 합은 길이 $P \/ 2$인 DFT임.
반으로 나눈 결과를 재사용하고 두 결과를 결합하므로 전체를 매번 새로 더할 필요가 없음.
이를 반복하는 것이 FFT의 기본 원리임.

== 12. 보충: 알려진 신호의 CWT 응답
<알려진-신호로-정답-유도하기>
이 절에서는 결과를 식으로 계산할 수 있는 신호를 사용해 CWT의 성질을 확인함.
DNA 신호와 별개인 수학적 예제도 포함하며, 실제 DNA 인코딩과의 차이는 12.4절에서 다룸.

=== 12.1 상수 신호의 응답
<상수-신호의-계수는-0입니다>
무한히 이어진 상수 $x \( t \) = c$에 대해서는

$ W_x \( a \, b \) = c sqrt(a) integral_(bb(R)) overline(psi \( u \)) thin d u = 0 . $

즉, 평균 0 조건 때문에 모든 스케일과 위치에서 계수가 0임.
그러나 유한 contig 전체가 A인 경우는 무한 상수와 다름.
contig 밖은 0이므로 양 끝에는 값의 변화가 있고 그 주변의 계수가 0일 필요는 없음.
절단된 커널이 contig 안에 전부 들어가는 내부 위치에서만 거의 0인지 확인함.
정확히 0이 아니라 거의 0인 이유는 꼬리 절단과 부동소수점 계산 오차가 있기 때문임.

=== 12.2 순수 정현파의 응답
<순수-정현파의-응답>
계산하기 쉬운 연속시간 신호 $x \( t \) = e^(i omega t)$를 넣음.
$t = b + a u$를 대입하면

$ W_x \( a \, b \) = sqrt(a) thin e^(i omega b) integral_(bb(R)) e^(i a omega u) overline(psi \( u \)) thin d u . $

$overline(psi \( u \)) = N e^(- u^2 \/ 2) \( e^(- 6 i u) - e^(- 18) \)$를 대입하면,
가우스 함수에 곱하는 진동의 각주파수는 각각 $a omega - 6$과 $a omega$가 됨.
11절의 가우스 적분을 두 번 사용하면

$
  #box(stroke: black, inset: 3pt, [$ W_x \( a \, b \) = N sqrt(2 pi a) thin e^(i omega b) [e^(- \( a omega - 6 \)^2 \/ 2) - e^(- 18) e^(- \( a omega \)^2 \/ 2)] . $])
$

이 식은 크기뿐 아니라 $e^(i omega b)$를 통해 복소 위상까지 알려줌.
켤레 방향이 틀리면 $a omega - 6$ 대신 다른 부호가 나타나므로 구현 검산에 유용함.
고정 스케일에서 주된 응답은 $a omega approx 6$인 입력 주파수에 나타남.
보정항이 매우 작으므로 중심 주파수와 대표 주기의 근삿값은

$
  omega_c approx 6 / a \, #h(2em) f_c approx frac(6, 2 pi a) \, #h(2em) T approx frac(2 pi a, 6) upright(" bp")
$

임. $T = 1 \/ f_c$이며 “스케일 $a$는 주기 $a$ bp”라고 해석하면 안 됨.
e.g.) 주기 3 bp에 맞는 스케일은 약 $a = 9 \/ pi approx 2.86$임. \
자연수 스케일 3은 가까운 패턴을 볼 수 있지만 정확히 같은 중심 주기는 아님.

c.f.) 3.1절의 시간-주파수 해상도 설명은 #text(weight: "bold")[동일한 분석 주파수에서 비교]할 때의 관계임. \
일반적인 $omega_0$를 사용하면 중심 주파수는 약 $omega_0 \/ a$임.
이를 고정한 채 $omega_0$를 키우면 $a$도 커져 창이 넓어지므로,
주파수 구분은 더 세밀해지고 위치 구분은 덜 정밀해짐.
반면 같은 $a$에서 $omega_0$만 바꾸면 가우스 창의 폭은 그대로이고 중심 주파수가 이동함.

스케일 1에서는 $f_c approx 0.955$ cycles/bp로,
간격 1인 샘플의 Nyquist 기준 $0.5$ cycles/bp보다 큼.

=== 12.3 작은 스케일이 유효하다는 말의 범위
<작은-스케일이-유효하다는-말의-범위>
정수 샘플에서는

$ e^(i \( omega + 2 pi k \) n) = e^(i omega n) quad \( k upright("는 정수") \) $

임. 서로 다른 연속 주파수가 같은 샘플을 만들 수 있다는 뜻임.
이 현상을 에일리어싱(aliasing)이라고 부름.
따라서 샘플만으로 원래 연속 주파수를 유일하게 알아낼 수 없음.
일반적인 샘플링 정리의 복원 보장에는 원 신호가 Nyquist 한계 안에 있다는 조건도 필요함.

현재 구현은 원 신호를 추측하는 대신 2절의 계단 신호를 선택해 그 CWT를 계산함.
계단 신호는 경계에서 불연속이므로 높은 주파수 성분도 포함함.
따라서 자연수 스케일 1의 적분도 정의됨.
이는 샘플에서 소실된 정보를 되찾았다는 뜻이 아니라,
#text(weight: "bold")[명시적으로 선택한 계단 신호를 분석함]이라는 뜻임.

=== 12.4 코드가 받는 정현파 샘플의 응답
<코드가-받는-정현파-샘플의-응답>
순수 연속 정현파와 그 샘플을 계단으로 재구성한 신호는 같지 않음.
`ACTG`를 반복하면 $x \[ n \] = e^(i pi n \/ 2)$이지만 구간 안에서는 정현파가 아니라 상수임.
따라서 실제 구현을 순수 정현파 공식과 그대로 비교하면 비교 기준 자체가 잘못됨.

이 계단 신호는 다음 회전 패턴들의 합으로 표현할 수 있음.
$omega = pi \/ 2$, $Omega_k = omega + 2 pi k$로 놓으면

$ x \( t \) = sum_(k in bb(Z)) c_k e^(i Omega_k t) \, #h(2em) c_k = frac(sin \( Omega_k \/ 2 \), Omega_k \/ 2) . $

이것은 주기 함수를 여러 주파수의 진동으로 분해하는 푸리에 급수(Fourier series)임.
$bb(Z)$는 정수 전체를 뜻하며 $k$가 변할 때 샘플에서 구분되지 않는 주파수들이 나타남.
각 성분의 크기는 폭 1인 구간을 유지한 재구성 규칙에 의해 $c_k$로 결정됨.
구간 내부와 적분에서 사용하는 전개이며, 점프 지점에서의 값은 계단 함수의 한쪽 값과 다를 수 있음.
그러나 점 하나의 값은 CWT 적분에 영향을 주지 않음.

계수는 각 폭 1 구간에서 다음 적분을 계산해 얻음.

$ integral_(- 1 \/ 2)^(1 \/ 2) e^(- i Omega_k t) thin d t = frac(2 sin \( Omega_k \/ 2 \), Omega_k) = c_k . $

길이 4인 한 주기에서 급수 계수는
$1 / 4 integral_(- 1 \/ 2)^(7 \/ 2) x \( t \) e^(- i Omega_k t) d t$임.
네 구간으로 나누고 $t = n + v$로 치환하면 $e^(i omega n) e^(- i Omega_k n) = 1$이므로
위 구간 적분이 네 번 더해지고 앞의 $1 \/ 4$와 상쇄됨.
그 외 길이 4의 푸리에 주파수에서는 네 복소수의 등비합이 0이 되어 계수도 0임.

CWT는 신호에 대해 선형이므로 각 회전 패턴의 응답을 더하면 됨.
선형이라는 것은 신호를 더한 뒤 변환한 결과와 각각 변환해 더한 결과가 같다는 뜻임.
정수 위치 $b$에서는 $e^(i Omega_k b) = e^(i omega b)$이므로

$
  #box(stroke: black, inset: 3pt, [$ W_x \( a \, b \) = e^(i omega b) N sqrt(2 pi a) sum_(k in bb(Z)) c_k [e^(- \( a Omega_k - 6 \)^2 \/ 2) - e^(- 18) e^(- \( a Omega_k \)^2 \/ 2)] . $])
$

가우스 항 때문에 중심에서 멀리 떨어진 주파수의 기여는 빠르게 작아짐.
따라서 충분한 범위의 유한 합으로 수치값을 근사해 실제 C 출력과 비교할 수 있음.
유한 contig에서는 커널이 모두 들어가는 내부 위치에서 비교해야 경계 효과를 제외할 수 있음.

=== 12.5 Wavelet의 적합성 조건
<웨이블릿-적합성-조건을-확인합니다>
푸리에 변환을
$hat(psi) \( xi \) = integral psi \( u \) e^(- i xi u) d u$로 정의하면

$ hat(psi) \( xi \) = N sqrt(2 pi) [e^(- \( xi - 6 \)^2 \/ 2) - e^(- 18) e^(- xi^2 \/ 2)] . $

$xi$는 Mother wavelet 좌표에서의 각주파수임.
$xi = 0$에서 두 항이 같아 $hat(psi) \( 0 \) = 0$임.
이는 평균 0 조건을 주파수 영역에서 다시 표현한 것임.
또한 이 함수는 매끄러우므로 0 근처에서 $\| hat(psi) \( xi \) \|$는 상수 곱하기 $\| xi \|$ 이하로 제한됨.
멀리서는 가우스처럼 빠르게 감소함. 따라서 표준 적합성 조건인

$ 0 < integral_(bb(R)) frac(\| hat(psi) \( xi \) \|^2, \| xi \|) thin d xi < oo $

가 성립함.
- 0 근처: 분자의 크기 제곱은 $\| xi \|^2$에 비례하는 상한을 가지므로,
  $\| xi \|$로 나눠도 적분항은 발산하지 않음.
- 무한대: 가우스 항의 빠른 감소 때문에 적분이 유한함.

이 조건은 연속 wavelet 분석의 이론적 조건임.
선택된 출력 스케일만으로 모든 입력을 완전히 복원할 수 있다는 보장은 아님.
역변환에는 신호의 범위, 분석할 스케일과 주파수 방향, 역변환 규칙을 별도로 정해야 함.

== 13. 보충: 수치 적분과 꼬리 오차
<gauss-legendre-계수-구성하기>
이 장에서는 8.3절의 수치 적분 위치·가중치를 구성하고, 8.2절의 꼬리 오차 상한을 유도함.
기본 미분과 다항식 계산을 사용하지만 고등학교의 통상적인 계산보다는 길어짐.

=== 13.1 다항식 적분을 맞추는 조건
<다항식-적분을-정확히-맞춥니다>
$\[ - 1 \, 1 \]$에서 적분을 $sum_j w_j f \( z_j \)$로 근사한다고 하자.
먼저 상수와 다항식의 적분을 정확히 맞추도록 조건을 설정할 수 있음.

$
  sum_j w_j z_j^m = integral_(- 1)^1 z^m thin d z = cases(delim: "{", 2 \/ \( m + 1 \) \, & m upright("은 짝수") \,, 0 \, & m upright("은 홀수") .)
$

e.g.) $m = 0$이면 $sum_j w_j = 2$, $m = 1$이면 $sum_j w_j z_j = 0$임. \
각각 상수 함수와 일차함수의 적분을 맞추는 조건임.

먼저 평가 위치가 2개인 경우를 생각하면 조건의 의미를 확인하기 쉬움.
두 위치를 $- s$, $s$로, 같은 가중치를 $w$로 놓음.
상수 함수와 이차함수의 적분을 맞추려면

$ 2 w = 2 \, #h(2em) 2 w s^2 = 2 / 3 $

이어야 함. 따라서 $w = 1$, $s = 1 \/ sqrt(3)$임.
일차·삼차함수는 대칭인 두 위치에서 값이 상쇄되므로 적분값 0과도 일치함.
결국 이 2점 구적법은 3차 이하의 다항식에 정확함.

$ integral_(- 1)^1 f \( z \) d z approx f \( - 1 \/ sqrt(3) \) + f \( 1 \/ sqrt(3) \) . $

16점 방법도 평가 위치와 가중치로 다항식의 적분을 맞춘다는 같은 원리임.
위치 16개가 이미 정해졌다면 $m = 0 \, dots.h \, 15$의 식은 가중치 16개에 대한 연립일차방정식임.
Gauss-Legendre 방법은 위치도 특별하게 선택해 정확히 맞는 차수를 31까지 높임.

=== 13.2 Legendre 다항식의 근으로 위치 선택
<legendre-다항식의-근을-위치로-사용합니다>
필요한 다항식은 다음 점화식으로 만들 수 있음.

$
  P_0 \( z \) = 1 \, quad P_1 \( z \) = z \, quad P_(m + 1) \( z \) = frac(\( 2 m + 1 \) z P_m \( z \) - m P_(m - 1) \( z \), m + 1) .
$

e.g.) $P_2 \( z \) = \( 3 z^2 - 1 \) \/ 2$임. 그 근은 앞의 2점 예제에서 구한 $plus.minus 1 \/ sqrt(3)$과 같음. \
$P_16 \( z \) = 0$의 16개 실근을 $z_j$로 선택함.
근을 수치적으로 구할 때는 부호가 달라지는 구간을 찾아 이분법을 사용하거나,
도함수를 이용한 Newton 방법을 사용할 수 있음.

같은 다항식은 다음 미분식으로도 표현됨.

$ P_m \( z \) = frac(1, 2^m m !) frac(d^m, d z^m) \( z^2 - 1 \)^m . $

이 식에서 16번 부분적분하면 차수가 15 이하인 다항식 $q$에 대해

$ integral_(- 1)^1 P_16 \( z \) q \( z \) thin d z = 0 $

을 얻음.
경계에서는 $\( z^2 - 1 \)^16$의 15차 이하 도함수가 0이라 경계항이 사라짐.
마지막에는 $q$의 16차 도함수가 0이 되므로 적분값도 0임.
이처럼 낮은 차수의 다항식과 곱해 적분하면 0이 되는 성질을 직교성(orthogonality)이라고 부름.
점화식과 미분식의 다항식이 같다는 것은 계수를 전개해 확인할 수 있음.

=== 13.3 31차까지 정확한 이유
<왜-31차까지-정확한가요>
차수가 31 이하인 다항식 $p$를 $P_16$으로 나누면

$ p \( z \) = q \( z \) P_16 \( z \) + r \( z \) \, #h(2em) deg q lt.eq 15 \, quad deg r lt.eq 15 . $

여기서 $deg$는 다항식의 차수를 뜻함.
- 적분: 앞 절의 직교성 때문에 $q P_16$ 부분은 0임.
- 구적합: $P_16 \( z_j \) = 0$이므로 역시 $q P_16$ 부분은 0임.
- 나머지 $r$: 15차 이하이며 가중치가 그 범위의 적분을 맞추도록 정해져 있음.

따라서 구적합과 적분이 같음. 이것이 31차까지 정확한 이유임.
단, 이는 정확한 위치와 가중치를 사용한 수학적 성질임.
코드에 저장한 소수 근삿값과 부동소수점 연산에는 작은 반올림 오차가 존재함.

실제 가중치는 다음 식으로도 구할 수 있음.

$ w_j = frac(2, \( 1 - z_j^2 \) \[ P'_16 \( z_j \) \]^2) . $

이 식 대신 13.1절의 연립방정식을 풀어도 같은 가중치를 구성할 수 있음.
대칭인 근에는 같은 가중치가 있으므로 코드에는 양의 근 8개만 저장함.
실수부와 허수부를 각각 같은 방법으로 적분하면 복소 적분을 계산할 수 있음.

=== 13.4 꼬리 오차 상한의 유도
8.2절에서는 꼬리를 버린 영향이 매우 작다는 결론을 사용했음.
그 근거는 wavelet 크기를 가우스 포락선으로 제한하는 것임.
삼각부등식으로 $\| e^(6 i u) - e^(- 18) \| lt.eq 1 + e^(- 18)$이므로

$ \| psi \( u \) \| lt.eq N \( 1 + e^(- 18) \) e^(- u^2 \/ 2) . $

$\| x \( t \) \| lt.eq 1$이고 양쪽 꼬리가 대칭이므로

$
  E_(upright(t a i l)) & lt.eq sqrt(a) integral_(\| u \| > 8) \| psi \( u \) \| thin d u \
                       & lt.eq 2 N \( 1 + e^(- 18) \) sqrt(a) integral_8^oo e^(- u^2 \/ 2) thin d u .
$

$u gt.eq 8$에서 $1 lt.eq u \/ 8$을 이용하면

$ integral_8^oo e^(- u^2 \/ 2) thin d u lt.eq 1 / 8 integral_8^oo u e^(- u^2 \/ 2) thin d u = e^(- 32) / 8 $

임. 이를 앞 식에 넣으면 8.2절의 상한을 얻음.
꼬리 적분을 특수함수로 쓸 수도 있음.

$ "erfc" \( v \) = 2 / sqrt(pi) integral_v^oo e^(- s^2) thin d s . $

$u = sqrt(2) s$로 치환하면 같은 상한을 다음과 같이 표현할 수 있음.

$ E_(upright(t a i l)) lt.eq 2 N \( 1 + e^(- 18) \) sqrt(a) sqrt(pi \/ 2) thin "erfc" \( 8 \/ sqrt(2) \) . $

== 14. 보충: C 구현의 공통 정의와 입출력
<구현과-수식의-대응>
앞 절에 임베드한 코드와 수식의 대응은 다음과 같음.
코드는 현재 `anno_cwt.c`의 함수 본문을 포함한 것으로, 단독 예제가 아니라 아래 공통 정의를 공유함.

#figure(
  align(center)[#table(
    columns: (50%, 50%),
    align: (auto, auto),
    table.header([구현 요소], [수학적 역할]),
    table.hline(),
    [`base_signal`], [염기를 $x \[ n \]$으로 변환],
    [`morlet_integral`],
    [$overline(psi \( u \))$의 구간 적분, 꼬리
      절단, 수치 구적],
    [`wavelets_init`],
    [각 스케일의 $k_a \[ j \]$, 반경과 커널 길이
      계산],
    [`convolution_fft_size`], [선형 합성곱에 충분한 FFT 길이 선택],
    [`cwt_extract`],
    [역순 커널과 입력의 FFT 곱, 역변환, 위치 이동,
      $1 \/ P$ 적용],
    [`export_cwt`], [두 가닥의 정수 위치 계수를 complex64로 저장],
  )],
  kind: table,
)

=== 14.1 공통 정의와 메모리 할당
`WAVE_COUNT`는 스케일 개수, `CWT_CHANNELS`는 실수부·허수부를 따로 센 채널 개수임.
`Wavelets`는 스케일별 커널과 길이를 저장함.
`MAX_WAVE_SIZE`는 초기 메모리 할당 크기의 기준이며 커널 절단 상한이 아님.
저장 공간은 각 커널에 맞춰 자동으로 확장되며, 초기 할당 기준이 0이나 음수이면 1을 사용함.
`alloc`은 `calloc`을 사용하므로 새 배열을 0으로 초기화하며, 크기 초과와 할당 실패를 확인함.
아래 헤더 및 자료형 정의가 앞의 함수들보다 먼저 필요함.

```c
#define _POSIX_C_SOURCE 200809L
#include "kseq.h"
#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

#define RFFT_IMPLEMENTATION
#include "rfft.h"

KSEQ_INIT(gzFile, gzread)

#ifndef CWT_SCALES
#define CWT_SCALES 4, 5, 6, 7, 8, 9
#endif
#ifndef MAX_WAVE_SIZE
#define MAX_WAVE_SIZE 9
#endif
static const double wave_scales[] = {CWT_SCALES};
#define WAVE_COUNT (sizeof(wave_scales) / sizeof(wave_scales[0]))
#define CWT_CHANNELS (2 * WAVE_COUNT)

typedef struct {
  double complex *kernel[WAVE_COUNT];
  int widths[WAVE_COUNT], max_width;
} Wavelets;
typedef struct { char *name, *seq; int length; } Contig;

#define fail(...) do { \
  fputs("Error: ", stderr); \
  fprintf(stderr, __VA_ARGS__); \
  fputc('\n', stderr); \
  exit(EXIT_FAILURE); \
} while (0)

static void *alloc(size_t count, size_t size) {
  if (size && count > SIZE_MAX / size) fail("Allocation size overflow");
  void *memory = calloc(count ? count : 1, size);
  if (!memory) fail("Out of memory");
  return memory;
}
```

=== 14.2 FASTA 입력과 바이너리 출력
FASTA는 `kseq.h`로 파싱하며 zlib을 통해 일반 파일과 gzip 압축 파일을 읽음.
`-`를 입력 경로로 주면 표준 입력을 사용함.
이 함수는 FASTA의 모든 contig를 메모리에 올림.
CWT 계산은 청크 단위이지만 입력 서열 전체를 스트리밍하는 구현은 아님.

```c
static Contig *read_fasta(const char *path, int *count) {
  gzFile input = gzopen(!strcmp(path, "-") ? "/dev/stdin" : path, "rb");
  if (!input) fail("Cannot open FASTA: %s", path);
  kseq_t *record = kseq_init(input);
  Contig *contigs = NULL;
  *count = 0;
  while (kseq_read(record) >= 0) {
    if (!record->seq.l) continue;
    if (record->seq.l > INT_MAX - 1024) fail("Contig too long: %s", record->name.s);
    Contig *grown = realloc(contigs, (size_t)(*count + 1) * sizeof(*contigs));
    if (!grown) fail("Out of memory");
    contigs = grown;
    Contig *contig = contigs + (*count)++;
    contig->name = strdup(record->name.s);
    contig->seq = strdup(record->seq.s);
    contig->length = (int)record->seq.l;
    if (!contig->name || !contig->seq) fail("Out of memory");
  }
  kseq_destroy(record);
  if (gzclose(input) != Z_OK || !*count) fail("No readable FASTA records: %s", path);
  return contigs;
}

static FILE *create_in(const char *directory, const char *name) {
  char *path = alloc(strlen(directory) + strlen(name) + 2, 1);
  sprintf(path, "%s/%s", directory, name);
  FILE *file = fopen(path, "wb");
  if (!file) fail("Cannot create %s", path);
  free(path);
  return file;
}

static void write_floats(FILE *stream, const float *values, size_t count) {
  _Static_assert(sizeof(float) == 4 && FLT_RADIX == 2 && FLT_MANT_DIG == 24,
                 "CWT export requires binary32 floats");
  uint16_t endian = 1;
  if (*(unsigned char *)&endian) {
    if (fwrite(values, sizeof(float), count, stream) != count) fail("Cannot write CWT matrix");
    return;
  }
  for (size_t index = 0; index < count; index++) {
    unsigned char bytes[4], swapped[4];
    memcpy(bytes, values + index, 4);
    for (int byte = 0; byte < 4; byte++) swapped[byte] = bytes[3 - byte];
    if (fwrite(swapped, 1, 4, stream) != 4) fail("Cannot write CWT matrix");
  }
}
```

`write_floats`는 실수부·허수부 순서를 유지하고,
실행 환경의 바이트 순서와 관계없이 little-endian으로 기록함.
즉, 파일을 읽는 쪽에서 동일한 complex64 형식으로 해석할 수 있도록 출력 규칙을 고정한 것임.

=== 14.3 실행 순서와 함수 연결
실행 순서는 다음과 같음.
+ `wavelets_init`: 스케일별 구간 적분 커널 생성.
+ `export_cwt`: FASTA 입력, 역상보 서열 생성, 양쪽 가닥을 청크 단위로 처리.
+ `cwt_extract`: 주변 문맥을 포함한 FFT 합성곱으로 계수 계산.
+ `write_floats`: complex64 형식으로 변환된 값을 파일에 기록.
+ 계산이 끝나면 입력 배열과 커널 메모리 해제.

```c
int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <genome.fa[.gz]> <new_directory>\n", argv[0]);
    return EXIT_FAILURE;
  }
  Wavelets wavelets;
  wavelets_init(&wavelets);
  export_cwt(&wavelets, argv[1], argv[2]);
  for (size_t scale = 0; scale < WAVE_COUNT; scale++) free(wavelets.kernel[scale]);
  return EXIT_SUCCESS;
}
```

문서의 코드 블록을 소스 파일로 연결할 때는 다음 순서가 필요함.
공통 정의 → 6.5절의 인코딩 → 8.3절의 적분 → 7.4절의 커널 생성 →
9.4절의 FFT 길이 선택 → 9.6절의 CWT 계산 → 14.2절의 입출력 → 10.4절의 저장 → `main`.
외부 헤더 `kseq.h`, `rfft.h`와 zlib은 별도로 필요함.


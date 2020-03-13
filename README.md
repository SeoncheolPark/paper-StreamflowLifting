# StreamflowLifting

- `code`: R code, directory path 바꿔야, 보안 문제로 구글 API 키 삭제.

- `data`: Geum-River dataset, shape 파일 (`River` 폴더 안에 있음) (혹시 모를 저작권 문제를 고려해 삭제)

	+ Data are gathered from the [Water Environment Information System](http://211.114.21.35/KRF_DEV/), operated by Ministry of Environment, South Korea.

	+ Geum-River network can be obtained in [Korean Reach File](http://water.nier.go.kr/front/riverNetwork/riverNetwork.jsp) in the Wter Information System, operated by National Institute of Environmental Research, South Korea.

- `result_RDS`: 시뮬레이션 result, 100번 iteration한 것을 `.RDS` 파일로 저장

- `stpca_package`: Flow-directed PCA 논문의 R 패키지를 그대로 가져왔다 (시뮬레이션 데이터 분석에 쓰임)

- `stream_result`: 논문에 쓰인 그림 저장

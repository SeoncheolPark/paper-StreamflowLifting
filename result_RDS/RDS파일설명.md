# 설명

- `GeumRiver(modified100trajectories).RDS`등: real data analysis nondecimated lifting scheme 100번 돌림, `GeumRiver.RDS`는 오리지널 버전, `GeumRiver(modified).RDS`는 장소 index를 바꿨을 때의 버전

- `ListStreamSim40(sd1)nlt_withresidual.RDS`등: simulation study의 미호천

- `ListStreamSTPCA80(sd1)nlt.RDS`등: simulation study의 `Gallacher`

- `RealListStreamSTPCA30(sd1)nlt_withresidual.RDS`등: simulationa study의 `Gallacher60`

- `ListStreamSim40(sd1)nlt_withresidualwithtimecorr.RDS`, `RealListStreamSTPCA60(sd1)nlt_withresidualwithtimecorr.RDS`등: simulation study에서 correlated error 생성시 결과물(revision 2, 이 때 random noise를 mean만 직전 장소에 dependent한 univariate normal로 만든 것은 defunct 버전으로 하고 conditional biavariate normal로 만드는 것으로 대체함)

- `weight_list.RDS`: real data analysis에서 `weight_vec_candidate`와 `weight_vec_candidate2`를 따로 저장해 두기 위해 생성함
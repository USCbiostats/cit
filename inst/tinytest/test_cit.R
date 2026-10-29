# Loading the test data
load("testpb.Rdata")

# Multiple runs give the same results
set.seed(1231)
ans1 <- cit.bp(L, G, T1, C)
ans2 <- cit.bp(L, G, T1, C)
ans3 <- cit.bp(L, G, T1, C)
ans4 <- cit.bp(L, G, T1, C)

expect_equal(ans1[2:4], ans2[2:4])
expect_equal(ans1[2:4], ans3[2:4])
expect_equal(ans1[2:4], ans4[2:4])

# The stochastic pvals should be roughly the same
expect_true(all(abs(ans1[c(1,5)] - ans2[c(1,5)]) < .04))
expect_true(all(abs(ans1[c(1,5)] - ans3[c(1,5)]) < .04))
expect_true(all(abs(ans1[c(1,5)] - ans4[c(1,5)]) < .04))


# Parametric test 4: a null statistic of exactly 0 with a central fit (cdf 0, z = -Inf)
# no longer makes the p-value NaN, and an observed statistic of 0 gives p > 0
t4p <- cit:::t4_para_p
st0 <- c(0, seq(0.001, 0.1, length.out = 49))
p0 <- t4p(0.002, st0, c(df1 = 1, df2 = 997))
expect_true(is.finite(p0) && p0 > 0 && p0 < 1)
expect_true(is.finite(t4p(0.002, st0, c(df1 = 1, df2 = NA))))
expect_true(t4p(0, st0, c(df1 = 1, df2 = 997)) > 0)

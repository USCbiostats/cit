
# Placeholder with simple test
expect_equal(1 + 1, 2)

load("testpb.Rdata")

# Multiple runs give the same results
ans1 <- cit.bp(L, G, T1, C)
ans2 <- cit.bp(L, G, T1, C)
ans3 <- cit.bp(L, G, T1, C)
ans4 <- cit.bp(L, G, T1, C)

expect_equal(ans1[2:4], ans2[2:4])
expect_equal(ans1[2:4], ans3[2:4])
expect_equal(ans1[2:4], ans4[2:4])

# The stochastic pvals should be roughly the same
set.seed(44)
expect_true(all(abs(ans1[c(1,5)] - ans2[c(1,5)]) < .02))
expect_true(all(abs(ans1[c(1,5)] - ans3[c(1,5)]) < .02))
expect_true(all(abs(ans1[c(1,5)] - ans4[c(1,5)]) < .02))


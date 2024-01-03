test_that("formula2terms", {
  f1 <- y ~ x1 + x2*x3
  f2 <- resid ~ s(checklist_number, observer_id, bs="fs")

  expect_setequal(formula2terms(f1), c('y', 'x1', 'x2', 'x3'))
  expect_setequal(formula2terms(f2), c('resid', 'checklist_number', 'observer_id'))
})


test_that('get_missing_dates errors on non-date class', {
  d1 <- c(1:5, 10)
  expect_error(get_missing_dates(d1))
})


test_that('names2chr',{
  expect_setequal(names2chr(foo, bar, buzz), c("foo", "bar", "buzz"))
})

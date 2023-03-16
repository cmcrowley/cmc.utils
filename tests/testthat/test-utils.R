test_that("formula2terms", {
  f1 <- y ~ x1 + x2*x3
  f2 <- resid ~ s(checklist_number, observer_id, bs="fs")

  expect_equal(formula2terms(f1), c('y', 'x1', 'x2', 'x3'))

  #TODO see mgcv::interpret.gam
  expect_equal(formula2terms(f2), c('checklist_number', 'observer_id'))
})


test_that('get_missing_dates errors on non-date class', {
  d1 <- c(1:5, 10)
  expect_error(get_missing_dates(d1))
})

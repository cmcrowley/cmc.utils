test_that("round_to_nearest", {

  # round to nearest base-x ---
  expect_equal(round_to_nearest(2.5, logbase=2),
               2)
  expect_equal(round_to_nearest(x = 2.5, logbase = 4),
               4)
  expect_equal(round_to_nearest(x=1, logbase=2),
               1)
  expect_equal(round_to_nearest(x=0, logbase=2),
               0)
  # Test x < 1
  expect_equal(round_to_nearest(x = 0.5, logbase=2),
               0.5)
  # Test x < 0
  expect_equal(round_to_nearest(x = -1, logbase=2),
               NA)


  # round to nearest decimal ---
  expect_equal(round_to_nearest(4.3, constant=0.5),
               4.5)

  # Test ties ---
  # (for now, not handled; always gets rounded up in log case and down in decimal case):
  expect_equal(round_to_nearest(96, logbase=2),
               128)
  expect_equal(round_to_nearest(2.25, 0.5),
               2)
})

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

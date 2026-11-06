
test_that("lod() catches bad `aff` input", {
  x = nuclearPed(1) |>
    addMarker()

  mod = diseaseModel("AD")

  expect_error(lod(x, aff = 1:4, mod = mod), "Invalid `aff`")
  expect_error(lod(x, aff = 0:1, mod = mod), "Unknown ID label: 0")
  expect_error(lod(x, aff = list(aff = 1), mod = mod), "Invalid name of `aff`")
  expect_error(lod(x, aff = list(affected = 1, unaffected = 1), mod = mod),
               "Individual with multiple affection statuses")
})


test_that("lod() with liability classes", {
  x = nuclearPed(3) |>
    addMarker(geno = c("1/2", "1/1", "1/2","1/2","1/2"))

  mod = diseaseModel("AD", pen = data.frame(f0 = 0, f1 = 1:0, f2 = 1:0))

  # Switch off last sib
  res = lod(x, aff = c(2,1,2,2,1), mod = mod, liab = c(1,1,1,1,2))

  expect_equal(round(res - log10(2), 3), 0)

  # Same: Set last sib to unknown
  expect_equal(res, lod(x, aff = c(2,1,2,2,0), mod = diseaseModel("AD")))
})

test_that("LOD is invariant to loop breaker choice", {
  x = cousinPed(1, child = TRUE) |>
    addMarker("9" = "2/2", alleles = 1:2, afreq = c(0.99, 0.01))
  mod = diseaseModel("AR", dfreq = 0.01)

  a = lod(x, aff = 9, model = mod, rho = 0.1, loopBreakers = 7)
  b = lod(x, aff = 9, model = mod, rho = 0.1, loopBreakers = 8)
  expect_equal(a, b)
})

test_that("peakSummary() adds physical positions", {
  x = data.frame(CHROM = 1, MARKER = c("a","b","c"),
                 MB = 1:3, LOD = c(0,3,0))
  pm = data.frame(MARKER = c("a","b","c"),
                  CHROM = 1, BP = (1:3)*1e6)

  z = peakSummary(lodPeaks(x, threshold = 2), physmap = pm)
  expect_equal(z$FROM_BP, 1e6)
  expect_equal(z$TO_BP, 3e6)
})

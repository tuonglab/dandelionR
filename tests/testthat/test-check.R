library(rlang)
library(methods)
library(BiocGenerics)

# Test for .classCheck
test_that(".classCheck works correctly", {
    expect_error(.classCheck(NULL, "character"), NA)
    expect_error(.classCheck("test", "character"), NA)
    expect_error(
        .classCheck(123, "character"),
        "The '123' must be character, not numeric"
    )
    expect_error(.classCheck(TRUE, "logical"), NA)
    expect_error(
        .classCheck("test", "logical"),
        "The 'test' must be logical, not character"
    )
})

# Test for .typeCheck
test_that(".typeCheck works correctly", {
    expect_error(.typeCheck(NULL, "character"), NA)
    expect_error(.typeCheck("test", "character"), NA)
    expect_error(
        .typeCheck(123, "character"),
        "The '123' must be character, not double"
    )
    expect_error(.typeCheck(TRUE, "logical"), NA)
    expect_error(
        .typeCheck("test", "logical"),
        "The 'test' must be logical, not character"
    )
})

# Additional tests for edge cases
test_that(".classCheck handles edge cases", {
    expect_error(.classCheck(list(), "list"), NA)
    expect_error(.classCheck(data.frame(), "data.frame"), NA)
    expect_error(.classCheck(matrix(), "matrix"), NA)
})

test_that(".typeCheck handles edge cases", {
    expect_error(.typeCheck(list(), "list"), NA)
    expect_error(.typeCheck(data.frame(), "data.frame"), NA)
    expect_error(.typeCheck(matrix(), "matrix"), NA)
})

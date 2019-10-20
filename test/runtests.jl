using Test
using Random

using StatsModelsTools
using StatsModels

@testset "Test Polynomial Formulas" begin

    data = (
            x=float.(1:10),
            y=float.(1:10).^3,
           )

    f = @formula(y ~ poly(x, 3))
    sch = apply_schema(f, schema(data))
    @inferred modelcols(sch.rhs, data)

    @test poly(data.x, 3) == data.x.^3
    @test modelcols(sch.rhs, data) == [data.x data.x.^2 data.x.^3]

end

@testset "Eval RMS Formulas" begin

    @test_throws Exception rcs([1,2,3], 0)
    @test_throws Exception rcs([1,2,3], 1)

    data = (
            x=(float.(1:10)),
            y=float.(1:10),
            z=collect(1:10),
           )

    xs = collect(0:10)
    @test rms_generate_knots(float.(0:100), 3) == [10, 50, 90]
    @test rms_generate_knots((0:100), 3) == [10, 50, 90]
    @test rms_generate_knots(vcat(float.(0:10), [missing]), 3) == [1, 5, 9]
    @test rms_generate_knots(xs, 3) == [1, 5, 9]
    @test rms_generate_knots(reverse(xs) |> collect, 3) == [1, 5, 9]
    @test rms_generate_knots(shuffle(xs) |> collect, 3) == [1, 5, 9]

    @test rcs(xs, [1, 5, 9]) == [0.0  0.0;
                                                  1.0  0.0;
                                                  2.0  0.015625;
                                                  3.0  0.125;
                                                  4.0  0.421875;
                                                  5.0  1.0;
                                                  6.0  1.921875;
                                                  7.0  3.125;
                                                  8.0  4.515625;
                                                  9.0  6.0;
                                                  10.0 7.5]

    xs_with_missing = [xs; missing]
    @inferred rms_generate_knots(xs, 3)
    @inferred rms_generate_knots(xs_with_missing, 3)
    @inferred rcs(xs_with_missing, [1, 5, 9])
    @test rms_generate_knots(xs_with_missing, 3) == [1, 5, 9]

    knots = [1, 5, 9]
    result = rcs(xs, knots);
    result_with_missing = rcs(xs_with_missing, knots)
    @test all(result_with_missing .=== [result; [missing missing]])

    f = @formula(y ~ rcs(x, 3))
    sch = apply_schema(f, schema(data))
    @inferred modelcols(sch.rhs, data)
    @test coefnames(sch)[2] == ["rcs(x, 3)[1]", "rcs(x, 3)[2]"]

end


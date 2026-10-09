module FirstUseWindTests
using Test, StaticArrays
import GRAMSuite
const G = GRAMSuite

# A Julia stand-in for dynamically included native bindings. No native library,
# kernel or atmosphere data is loaded. Define the bindings after entering this
# function, then query in the same ordinary caller frame (no outer invokelatest).
function cold_query(; mode="perturbed", wind=true, present=true, nonfinite=false)
    binding = Module(gensym(:ColdWindBinding))
    trace = Any[]
    call_lock = ReentrantLock()
    winds = nonfinite ?
        (ewWind=NaN, nsWind=4.0, verticalWind=Inf,
         perturbedEWWind=30.0, perturbedNSWind=NaN, perturbedVerticalWind=50.0) :
        (ewWind=3.0, nsWind=4.0, verticalWind=5.0,
         perturbedEWWind=30.0, perturbedNSWind=40.0, perturbedVerticalWind=50.0)
    Core.eval(binding, quote
        struct Atmosphere end
        const TRACE = $trace
        const CALL_LOCK = $call_lock
        set_position!(a; kwargs...) = push!(TRACE,
            (event=:position, locked=islocked(CALL_LOCK), options=(; kwargs...)))
        update!(a) = (push!(TRACE, (event=:update, locked=islocked(CALL_LOCK))); 0)
        get_dynamics_state(a) = (density=1.25, temperature=211.0)
    end)
    if present
        Core.eval(binding, :(get_winds_state(a) = $winds))
    end
    old_frame_visible = isdefined(binding, :get_winds_state)
    latest_visible = Base.invokelatest(isdefined, binding, :get_winds_state)
    atmosphere = Base.invokelatest(Base.invokelatest(getfield, binding, :Atmosphere))
    model = G.GRAMAtmosphereModel(binding, atmosphere, "", "", "", "mars",
        G.InitialTime(), true, "")
    result = withenv("SPACEAGORA_GRAM_WIND_MODE" => mode) do
        G.point_density_state(model, 150_000.0, deg2rad(45.0), deg2rad(-20.0),
            7.0, wind; lock_obj=call_lock)
    end
    return (; result, trace, old_frame_visible, latest_visible, unlocked=!islocked(call_lock))
end

function check_query(result, expected)
    @test result.result == (1.25, 211.0, SVector{3,Float64}(expected))
    @test [event.event for event in result.trace] == [:position, :update]
    @test all(event.locked for event in result.trace)
    @test result.unlocked
    @test result.trace[1].options.height == 150.0
    @test result.trace[1].options.latitude ≈ 45.0
    @test result.trace[1].options.longitude ≈ -20.0
    @test result.trace[1].options.elapsed_time == 7.0
    @test result.latest_visible
    # Strict Julia 1.12 rejects stale global-binding access. Earlier supported
    # versions still exercise the ordinary-call API contract.
    if VERSION >= v"1.12" && Base.JLOptions().depwarn == 2
        @test !result.old_frame_visible
    end
end

const ORIGINAL_EPHEMERIS = G._GRAM_EPHEMERIS_STATE_FN[]
const ORIGINAL_MISSING_WARNING = G._GRAM_WIND_WARNING_EMITTED[]
const ORIGINAL_NONFINITE_WARNING = G._GRAM_NONFINITE_WIND_WARNING_EMITTED[]
try
    G._GRAM_EPHEMERIS_STATE_FN[] = nothing
    @testset "first-use wind binding visibility" begin
        for (mode, wind, expected) in (("perturbed", true, (30.0,40.0,50.0)),
                                      ("perturbed", false, (3.0,4.0,5.0)),
                                      ("nominal", true, (3.0,4.0,5.0)),
                                      ("nominal", false, (3.0,4.0,5.0)))
            G._GRAM_WIND_WARNING_EMITTED[] = false
            result = @test_logs cold_query(; mode, wind)
            check_query(result, expected)
            @test !G._GRAM_WIND_WARNING_EMITTED[]
        end
        @testset "existing nonfinite component policy" begin
            for (wind, expected) in ((true, (30.0,0.0,50.0)), (false, (0.0,4.0,0.0)))
                G._GRAM_NONFINITE_WIND_WARNING_EMITTED[] = false
                result = @test_logs (:warn, r"non-finite wind component") cold_query(; wind, nonfinite=true)
                check_query(result, expected)
                @test G._GRAM_NONFINITE_WIND_WARNING_EMITTED[]
            end
        end
        @testset "genuinely absent binding" begin
            G._GRAM_WIND_WARNING_EMITTED[] = false
            result = @test_logs (:warn, r"does not expose wind state") cold_query(; present=false)
            @test result.result == (1.25,211.0,SVector(0.0,0.0,0.0))
            @test !result.latest_visible
            @test [event.event for event in result.trace] == [:position, :update]
            @test all(event.locked for event in result.trace)
            @test result.unlocked
            @test G._GRAM_WIND_WARNING_EMITTED[]
            @test_logs cold_query(; present=false) # warning remains once-only
        end
    end
finally
    G._GRAM_EPHEMERIS_STATE_FN[] = ORIGINAL_EPHEMERIS
    G._GRAM_WIND_WARNING_EMITTED[] = ORIGINAL_MISSING_WARNING
    G._GRAM_NONFINITE_WIND_WARNING_EMITTED[] = ORIGINAL_NONFINITE_WARNING
end
end

using Test
# Pkg.test may use a compatibility mode that permits stale binding access.
# Also exercise the same strict mode as SpaceAGORA's loading diagnostics.
if Base.JLOptions().depwarn != 2
    @testset "strict first-use wind subprocess" begin
        command = `$(Base.julia_cmd()) --startup-file=no --history-file=no --depwarn=error --project=$(dirname(Base.active_project())) $(@__FILE__)`
        output = IOBuffer()
        process = run(pipeline(ignorestatus(command); stdout=output, stderr=output))
        text = String(take!(output))
        success(process) || print(stderr, text)
        @test success(process)
    end
end

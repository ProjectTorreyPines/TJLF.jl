using Test
using TJLF

directory = joinpath(@__DIR__, "test_eigen")

@testset "$directory" begin
  inputTJLF = readInput(joinpath(directory,"input.tjlf"))
  TJLF.run_tjlf(inputTJLF)
  fluxesJulia = sum(TJLF.run_tjlf(inputTJLF); dims=1)[1, :, :]

  # took these results from previous implementation with Arpark:eigs
  Fluxes_results=[-39.560547433921755, -34.280893518792716, -0.7300925138635125, 240.91550661719268, 180.22103212060173, 5.44996230222146, 0.010288609198657303, 5.661387914783233, 0.7838629459587043]

  for i in 1:3*inputTJLF.NS
    @test isapprox(sum(fluxesJulia[i]), sum(Fluxes_results[i]), rtol=5e-3)  
  end


  satParams = get_sat_params(inputTJLF)

  @test isapprox(1.0, satParams.SAT_geo0, rtol=1e-6)
  @test isapprox(0.4271492302106564, satParams.SAT_geo1, rtol=1e-6)
  @test isapprox(0.2726096564842753, satParams.SAT_geo2, rtol=1e-5)
  @test isapprox(3.855959015136675, satParams.R_unit, rtol=1e-6)
  @test isapprox(0.508492296295847, satParams.Bt0, rtol=1e-6)
  @test isapprox(1.15578209364075, satParams.grad_r0, rtol=1e-6)
end


#---------------------------------------------------------------------------------------
#------------------------------------------------------------------------
#------ Redo several regression tests with KrylovKit eigensolver, I randomly choose test01, test05, test06, test08, test10, test11
#----------------------------------------------------------------------
#--------------------------------------------------------------------------
# tglf regression test
directory = joinpath(@__DIR__, "tglf_regression")
tests = readdir(directory)

# 03 is s-alpha geometry, not implemented
excludeFolders = ["tglf02","tglf03","tglf04","tglf07","tglf09"]
testFolders = [joinpath(directory, item) for item in readdir(directory) if isdir(joinpath(directory, item)) && item ∉ excludeFolders]

for baseDirectory in testFolders
    @testset "$baseDirectory" begin
     
        fileDirectory = joinpath(baseDirectory, "out.tglf.gbflux")
        lines = readlines(fileDirectory)
        fluxesFortran = parse.(Float64, split(lines[1]))
        
        inputTJLF = readInput(joinpath(baseDirectory,"input.tglf"))
        inputTJLF.FIND_EIGEN =true
        inputTJLF.SMALL = 0.0
        fluxesJulia = sum(TJLF.run_tjlf(inputTJLF); dims=1)[1, :, :]
        
     
        for i in 1:3*inputTJLF.NS
            @test isapprox(sum(fluxesJulia[i]), sum(fluxesFortran[i]), rtol=1e-3, atol=1e-10)
        end
        for i in 3*inputTJLF.NS+1:4*inputTJLF.NS
            @test isapprox(sum(fluxesJulia[i+inputTJLF.NS]), sum(fluxesFortran[i]), rtol=1e-3, atol=1e-10)
        end
    end
end


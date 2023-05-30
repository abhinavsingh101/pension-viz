### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 620e83fe-0e29-4217-85eb-7cb7284ef499
using CSV, DataFrames, Plots, PlutoUI, PlotThemes

# ╔═╡ 64eb2201-2af5-4b06-bc76-cef60ba0cad4
begin
	df = CSV.read("mortality_data.csv", DataFrame);
	df_subset = select(df,[:Age,:Value,:Time]);
	replace!(df_subset.Age, "100+" => "100");
	df_subset.Age = map(x -> parse(Int, coalesce(x,"")), df_subset.Age);
	df_filtered = filter(row -> row.Time == 2019, df_subset);
	mortality_rate = select(df_filtered,[:Age,:Value]);
	rename!(mortality_rate, :Value => :probability_death);
end

# ╔═╡ 2fd5e8fe-c195-4a71-8f42-2260bad5ce0a
function m_rate_generalized(age::Int, Time, survival_table=mortality_rate)
	t = Symbol("Time_$Time")
	return survival_table[survival_table.Age .== age, t][1]
end

# ╔═╡ 9a00fe06-9828-4208-b655-58e933aebc08
function calculate_survivorfn_generalized(age::Int, Time, survival_table=mortality_rate)
    cooked = fill(NaN, 101)
    cooked[age+1] = 100
    for i in (age+2):100
        cooked[i] = cooked[i-1] - m_rate_generalized(i, Time, survival_table)*cooked[i-1]
    end
    return cooked
end

# ╔═╡ ab4bd33c-9e04-427d-8808-fdfc2ea56742
function m_rate(age::Int) 
	return mortality_rate[mortality_rate.Age .== age, :probability_death][1]
end

# ╔═╡ 3de95442-fca6-11ed-17dd-41fc8582cd5c
function calculate_survivorfn(age::Int)
    cooked = fill(NaN, 101)
    cooked[age+1] = 100
    for i in (age+2):100
        cooked[i] = cooked[i-1] - m_rate(i)*cooked[i-1]
    end
    return cooked
end

# ╔═╡ b27f4c1d-8ee4-4d0b-97ab-380d89dacd28
# begin
# 	mortality_rate.cooked_60 = calculate_survivorfn(60)
# 	mortality_rate.cooked_35 = calculate_survivorfn(35)
# end

# ╔═╡ 06384b92-5333-4b5f-849a-fb13ff6b571e
for a in 35:99
    col_name = Symbol("cooked_$a")
    mortality_rate[!, col_name] = calculate_survivorfn(a)
end

# ╔═╡ 2b76ad55-d09f-4fa8-9a98-ea3ff32b3ac6
# begin
# 	using Plots

# 	# Set the margins of the plot
# 	plot!(layout=grid(1, 1, heights=[0.8], widths=[0.8, 0.2]))
	
# 	# Plot the survivor function starting at age 35
# 	plot!(35:100, mortality_rate[36:101,:cooked_35], linewidth=2, color=:blue, xlabel="Age", ylabel="Survivors", label="Starting at age 35")
	
# 	# Plot the survivor function starting at age 60
# 	plot!(60:100, mortality_rate[61:101,:cooked_60], linewidth=2, color=:red, xlabel="Age", ylabel="Survivors", label="Starting at age 60")
	
# 	# Add a legend to the plot
# 	plot!(legend=:bottomleft, linewidth=2, color=[:blue, :red], framestyle=:none, label=["Starting at age 35" "Starting at age 60"])

# end

# ╔═╡ ebd0ce86-e18f-4d22-bcdf-d88c2b9d5b3f
function value_of_pension(p, S, r, l=0)
    sum([(1/r)^i * p[i] * S[i] for i in 1:length(p)][1:l])
end

# ╔═╡ 277797ab-724a-4848-afb3-1d2265ffde48
# function value_of_pension2(p, S, r, l=0)
#     n = min(l, size(p, 1))
#     PV = 0
#     for j in 1:size(p, 2)
#         PV += sum((1 ./ r) .^ (1:n) .* (p[1:n, j] .* S[1:n]))
#     end
#     return PV
# end

# ╔═╡ feb4581a-c97e-4259-a460-531185d28326
function value_of_pension2(p, S, r, l=0)
    n = min(l, size(p, 1))
    PVs = zeros(3)
    
    for j in 1:size(p, 2)
        PVs[j] = sum((1 ./ r) .^ (1:n) .* (p[1:n, j] .* S[1:n]))
    end
    
    return PVs
end

# ╔═╡ b331d2ba-67da-44cd-8101-f1ec6905d77f
# function value_of_pension(p, S, r, l=0)
#     discounted_value = [(1/r)^i * p[i] * S[i] for i in 1:length(p)]
#     if l > 0
#         discounted_value = discounted_value[1:l]
#     end
#     return sum(discounted_value)
# end

# ╔═╡ ae898e69-1ec1-4ca1-943c-99a8af3096cc
#As an example: Price an annuity of Rs.1 per day at age 60

# ╔═╡ f3579a67-2d9b-4e18-91b9-54aa806c4418
value_of_pension(fill(365, 40), mortality_rate[61:100,:cooked_60] ./ 100, 1.07, 40)

# ╔═╡ ea710959-34dd-49f4-851b-a029c045c757
function make_pension_mat(l, index,pension_amount)
    p1 = fill(NaN, (l, 3))
    p1[1,:] .= pension_amount
    for i in 2:l
        p1[i,1] = p1[i-1,1] + index[1] * p1[i-1,1]
        p1[i,2] = p1[i-1,2] + index[2] * p1[i-1,2]
        p1[i,3] = p1[i-1,3] + index[3] * p1[i-1,3]
    end
    return p1
end

# ╔═╡ 0f5333ef-3039-4470-8886-190ee368dd92


# ╔═╡ 0d141036-2508-4fc5-a0fd-f23109638f3f


# ╔═╡ f566dfb2-66da-40f3-bf03-0992d769873d


# ╔═╡ 2abd32aa-04b1-46e3-bf5b-740530773a2c
begin
	# Define the inflation and wage growth scenarios as arrays
	price_index = [0.03, 0.04, 0.05]
	wage_index = [0.07, 0.08, 0.09]
end

# ╔═╡ 60dbedcb-01f1-454b-9146-0a40785be747
if indexation == "prices"
    index = price_index
elseif indexation == "wages"
    index = wage_index
else
    error("Invalid indexation type")
end

# ╔═╡ 6914def4-a379-4b3b-8ebe-d307b7d3ca1f
# col = Symbol("cooked_$retirement_age")

# ╔═╡ 5370a609-07d7-4ef5-8dbe-361e25fbe853
# ╠═╡ disabled = true
#=╠═╡
begin
	# Price indexation
	p_mat = make_pension_mat(years_after_retirement, index)
	
	
end
  ╠═╡ =#

# ╔═╡ cf4e0c5d-2bdb-4a1c-895d-e7d3570627c7


# ╔═╡ 4030a33b-5771-4fb5-b8f8-53496cfe13ad
begin
	# Define the inflation and wage growth scenarios as arrays
	inflation = [0.03, 0.04, 0.05]
	wagegrowth = [0.07, 0.08, 0.09]
end

# ╔═╡ 06351924-c4d1-4ab8-bde3-3a3504677244
make_pension_mat(40, inflation, 365);

# ╔═╡ 4815ae30-f791-4bb3-8543-8d96193c461e
begin
	# Price indexation
	p1 = make_pension_mat(40, inflation, 365)
	value_of_pension2(p1, mortality_rate[61:100,:cooked_60] ./ 100, 1.07, 40)
	
end

# ╔═╡ cb87055a-a240-4129-98f6-6d87a6bdd8fa
begin
	# Wage indexation
	p2 = make_pension_mat(40, wagegrowth, 365)
	value_of_pension2(p2, mortality_rate[61:100,:cooked_60] ./ 100, 1.07, 40)
end

# ╔═╡ e7d550b2-be57-41d9-89aa-ffa4980eae3f
value_of_pension2(p1, mortality_rate[61:100,:cooked_60] ./ 100, 1.07, 40)

# ╔═╡ 345ce785-fc32-4340-8d42-e557bfe5d3f1
# begin	
# 	# Values for x (age) and y (survival probability)
# 	x = collect(60:100)
# 	y = mortality_rate[61:100,:cooked_60] ./ 100
	
# 	# Plot the survival curve
# 	plot(x, y, xlabel="Age", ylabel="Survival Probability", 
# 	     title="Survival Curve Starting at Age 60", linewidth=2)
	
# 	# Define scenarios for inflation and wage growth
# 	inflation_scenarios = [0.03, 0.04, 0.05] 
# 	wagegrowth_scenarios = [0.07, 0.08, 0.09]
	
# 	# Get NPV for inflation indexed annuity
# 	p_mat_inflation = make_pension_mat(40, inflation_scenarios, 365)
# 	npv_infl = value_of_pension2(p_mat_inflation, y, 1.07, 40)
	
# 	# Get NPV for wage indexed annuity
# 	p_mat_wage = make_pension_mat(40, wagegrowth_scenarios, 365) 
# 	npv_wage = value_of_pension2(p_mat_wage, y, 1.07, 40)
	
# 	# Plot the NPVs
# 	bar([npv_infl], label="Inflation Indexed")
# 	bar!([npv_wage], label="Wage Indexed")
# 	xlabel!("Scenario")
# 	ylabel!("NPV (Rs.)")
# 	title!("NPV of Annuities Starting at Age 60")
	
# 	# Add legend
# 	# legend!(:topright)
# end

# ╔═╡ c9351453-1d94-4360-b393-4e74b0640ba9



# ╔═╡ fb427276-0d24-4570-aced-f457c9477b02
# md"""
# $(@bind retirement_age Slider(35:60))
# $(@bind pension_amount Slider(365:1000))
# $(@bind c Button("🎖️"))
# $(@bind indexation Select(["prices" => "🍞", "wages" => "🧑🏽‍🤝‍🧑🏽"]))
# """

# ╔═╡ b4b4ca0d-8a1f-4ea2-8250-73015865c226
# print("""Lets calculate the present value of the pension for someone retiring at $retirement_age.
# Their pension will be Rs $pension_amount a year.
# The pension will be indexed to $indexation
# The maximum expected lifespan is $max_expected_lifespan""")

# ╔═╡ 82d54039-bca5-4c3b-8ef1-8ad88958c777
begin
	m1 = CSV.read("mortality_data_history.csv", DataFrame)
	m2 = select(m1,[:Age,:Value,:Time])
	replace!(m2.Age, "100+" => "100")
	m2.Age = map(x -> parse(Int, coalesce(x,"")), m2.Age)
	rename!(m2, :Value => :probability_death)
end

# ╔═╡ 4ae1dbcc-0938-4605-91e6-9c0da8a4c85b
# m2
reshaped_m2 = DataFrame(Age = unique(m2.Age));

# ╔═╡ 3704356c-44fd-4bd3-87f8-3652f83de810
unique_times = unique(m2.Time);

# ╔═╡ d5ec763f-6021-4566-a16f-4c219fe5116e
for time in unique_times
    time_col = m2[m2.Time .== time, :probability_death]
    time_col_name = Symbol("Time_", string(time))
    reshaped_m2[!, time_col_name] = time_col
end

# ╔═╡ f930becd-d683-434f-b97f-3fe81a1f8e8b
reshaped_m2

# ╔═╡ 185ae1c3-a68b-415d-b39c-be751c0779ff
for year in 1950:2019
    col_name = Symbol("cooked_60_$year")
    reshaped_m2[!, col_name] = calculate_survivorfn_generalized(60, year, reshaped_m2)
end

# ╔═╡ 92d52b8e-3436-45b6-95ac-6c57602d3c30
# begin	
# 	# Values for x (age) and y (survival probability)
# 	x = collect(61:100)
# 	y = reshaped_m2[61:100,col_viz] ./ 100
	
# 	# Plot the survival curve
# 	plot(x, y, xlabel="Age", ylabel="Survival Probability", 
# 	     title="Survival Curve Starting at Age 60 in 1950", linewidth=2)

# 	plot(x, y, xlabel="Age", ylabel="Survival Probability", 
# 	     title="Survival Curve Starting at Age 60 in 2019", linewidth=2)
# end

# ╔═╡ c40a18ce-5bfe-4ad0-91b0-b0454136abf3


# ╔═╡ d6344906-be05-4f49-9f19-9b770cd36d35


# ╔═╡ a077cc3b-d183-4528-b5ba-3c924bf929f1


# ╔═╡ babc50fd-c355-4ea9-a749-eb206d0e3495
# line_colors = cgrad(:Blues, 3, categorical=true)

# ╔═╡ ba0b0292-8188-4797-8185-95ffc29b0d25


# ╔═╡ a804f388-31de-43db-ba0e-600d02c26f06
function label_points(x_values, y_values, point_color, label_color, halign=:right, valign=:top,hpadding=10, vpadding=10)
    for (x_val, y_val) in zip(x_values, y_values)
        scatter!([x_val], [y_val], markershape=:circle, markercolor=point_color, markerstrokealpha = 0.2, markersize=6, label=nothing)
        annotate!([(x_val, y_val, Plots.text(string(round(y_val * 100, digits=1)) * "%", 8, halign, valign, color=label_color, hpadding, vpadding))])
    end
end

# ╔═╡ 73a8b5ea-853a-41db-85f8-57ab8936996c


# ╔═╡ c27ee5c1-997f-44a9-bc40-876bdc320e55


# ╔═╡ 0cb6f88d-7b8a-4d82-9795-310563aa31b4


# ╔═╡ 8a7bdd87-9a3a-404e-8a16-71ead8074072
# function cost_of_pension_over_time(;pension_amount=365, ylims=(1000, 35000), yticks=0:5000:25000)

# 	col_viz2 = Symbol("cooked_60_$year_viz")
# 	cost_var = value_of_pension(fill(pension_amount, 40), reshaped_m2[61:100,col_viz2] ./ 100, 1.07, 40)
# 	cost_1950 = value_of_pension(fill(pension_amount, 40), reshaped_m2[61:100,:cooked_60_1950] ./ 100, 1.07, 40)
# 	cost_2019 = value_of_pension(fill(pension_amount, 40), reshaped_m2[61:100,:cooked_60_2019] ./ 100, 1.07, 40)

# 	# Create a bar plot with the three pension costs
# 	fig = bar(["1950", "$year_viz", "2019"], [cost_1950, cost_var, cost_2019],
#           legend=false, ylims=(1000, 35000), yformatter=:plain, yticks=0:5000:25000,
#           color=[line_colors[1], line_colors[2], line_colors[3]], size=(300, 600))
	
# 	# Add plot title and axis labels
# 	title!("Pension Costs Comparison",titlefont=font(12), titlefontcolor=:white)
# 	xlabel!("Scenarios")
# 	ylabel!("Pension Costs")
# 	return fig
# end

# ╔═╡ a7a17973-c45b-4826-b156-b20f4ea28b2f


# ╔═╡ 6f752523-f2e2-4321-9453-3de46d5b4169
md"""
# Begins
"""

# ╔═╡ ae181d07-096c-4428-831a-5ee805810fef
md"""
We all know that life expectancy has increased over the years. This is one of the remarkable achievements of the 20th century.

But a longer retirement needs financial planning. If the gopvrnment **guranatees** a pension, then it certinaly needs to think about the cost of making that promise.

Life expectancy matters to pensions because when people live longer, their lifetime pensions are greater than if they die soon after retiring. Therefore, the pension burden of guaranteed pensions has increased over time.

By how much? We will find out.

But first, take a look at the survival probabilities of someone retiring at age 60."""

# ╔═╡ 15e4be92-93a7-4273-9e32-856006ddcb61
md"""

This graph shows the probability of survival in each year after retirement.

If you retire on your 61st birthday, the probability that you are alive is 100%. After that, it keeps decreasing every year.

Notice the numbers. They show the probability that you would live to a certain age if you retired at the age of 60 in 1950. The graph suggests a:

- 61.8% chance that you would have lived to age 70;
- 18.7% chance that you would live to 80; and
- 1% chance that you would have lived to 90

"""

# ╔═╡ 813f506d-1777-49a7-afd7-e962588b6290
md"""

Check this to add the survival rates for 2019. $(@bind enable_2019 CheckBox(default=false))

"""

# ╔═╡ 6be6c3c0-3974-42d6-a99e-213fd948bda5
md"""
Notice the differences.

For someone retiring at 60 in 2019, what is the probability that they would live to 90?
$(@bind question TextField(default="%"))
"""

# ╔═╡ bc867c54-b9b4-4b72-97c8-85af8136afb0
if question == "9.6%"
	println("Correct")
else
	println("It is 9.6%")
end

# ╔═╡ 5213433c-9913-4d83-8b9e-83c09fbac19c
md"""

Want to see survival rates for another year? $(@bind enable_var_year CheckBox(default=false))

"""

# ╔═╡ 51a13f1a-d305-4379-a23d-dc1affabeb92
if enable_var_year
    @bind year_viz Slider(1950:2019, default=1950, show_value=true)
else
    year_viz = 1950
    HTML("") # Display an empty HTML element
end

# ╔═╡ 222adafb-24bb-443c-be0c-bee6d69f641d
begin
	theme(:dark)
	# Values for x (age) and y (survival probability)
	col_viz = Symbol("cooked_60_$year_viz")
	x = collect(61:100)
	y = reshaped_m2[61:100,col_viz] ./ 100
	y_2019 = reshaped_m2[61:100,:cooked_60_2019] ./ 100
	y_1950 = reshaped_m2[61:100,:cooked_60_1950] ./ 100

	# Define the color palette
	point_colors = cgrad(:Blues, 3, categorical=true)
	line_colors = cgrad(:GnBu_3, 3, categorical=true)
	label_colors =cgrad(:Greens_3, 3, categorical=true) 
	
	# Plot the survival curve
	fig1 = plot(x, y_1950, linewidth=2, label = "1950",color=line_colors[1])

	# Plot the y-axis values for specific x-axis values for each line
	x_values = [70, 75, 80, 85, 90]  # Specify the x-axis values

	# Get the corresponding y-axis values for each line
	y_values_current = [y[ifelse(findfirst(isequal(x_val), x) === nothing, 1, findfirst(isequal(x_val), x))] for x_val in x_values]
	y_values_1950 = [y_1950[ifelse(findfirst(isequal(x_val), x) === nothing, 1, findfirst(isequal(x_val), x))] for x_val in x_values]
	y_values_2019 = [y_2019[ifelse(findfirst(isequal(x_val), x) === nothing, 1, findfirst(isequal(x_val), x))] for x_val in x_values]
	
	# Add labels to each point on each line

	label_points(x_values, y_values_1950, point_colors[2], label_colors[2],:right, :top)
	
	if enable_var_year
		plot!(fig1, x, y, linewidth=2, label = "$year_viz",color=line_colors[2], size=(800, 600))			
		label_points(x_values, y_values_current, point_colors[1], label_colors[1],:center, :top)

	end

	if enable_2019
		plot!(fig1, x, y_2019, linewidth=2, label = "2019",color=line_colors[3])
		label_points(x_values, y_values_2019, point_colors[3], label_colors[3],:left, :bottom)
	end
	
	
	plot!(xlabel="Age", ylabel="Survival Probability", title="Survival Curve Starting at Age 60", size=(800, 700))
end

# ╔═╡ 072efd53-a9e5-4ad5-bef3-bb3bcb33b3eb
line_colors

# ╔═╡ 4d22d4c3-ef7e-4e35-8f95-645830f8929b
point_colors

# ╔═╡ f41dcd9f-9611-4ed5-92ba-6c5009221fcb
fig1

# ╔═╡ 55da71d5-f553-4b05-9109-af83eadbaad6
cost_1950, cost_var, cost_2019, year_viz

# ╔═╡ e28f31bd-d229-4e0e-a48f-6854493fedf5
function cost_of_pension_over_time(;pension_amount=365, ylims=(1000, 35000), yticks=0:5000:25000)
    col_viz2 = Symbol("cooked_60_$year_viz")
    cost_var = value_of_pension(fill(pension_amount, 40), reshaped_m2[61:100,col_viz2] ./ 100, 1.07, 40)
    cost_1950 = value_of_pension(fill(pension_amount, 40), reshaped_m2[61:100,:cooked_60_1950] ./ 100, 1.07, 40)
    cost_2019 = value_of_pension(fill(pension_amount, 40), reshaped_m2[61:100,:cooked_60_2019] ./ 100, 1.07, 40)

    # Create a bar plot with the three pension costs
    fig = bar(["1950", "$year_viz", "2019"], [cost_1950, cost_var, cost_2019],
               legend=false, ylims=ylims, yformatter=:plain, yticks=yticks,
               color=[line_colors[1], line_colors[2], line_colors[3]], size=(300, 600))

    # Add plot title and axis labels
    title!("Pension Costs Comparison", titlefont=font(12), titlefontcolor=:white)
    xlabel!("Scenarios")
    ylabel!("Pension Costs")
    return fig
end


# ╔═╡ f040371c-f4a7-467c-90cb-d32ede5c6d2f
md"""
The survival rates are an important variable when thinking about pensions. When there is a greater chance that a pensioner will live to be 90 years old, the liability oif that pension increases.

Lets see how the cost of a modest pension of Rs 365 per year has changed over the years.
"""

# ╔═╡ 2667661c-4806-4a3a-8e56-3691d2ac93ce
begin
    pension_cost_comparison = cost_of_pension_over_time(pension_amount=365, ylims=(1000, 5000), yticks=0:500:5000)
end

# ╔═╡ 2b046c22-5b07-4b31-8ec5-351817273a0c
begin
	col_viz3 = Symbol("cooked_60_$year_viz")
	scenario_var = value_of_pension(fill(365, 40), reshaped_m2[61:100,col_viz3] ./ 100, 1.07, 40)
	scenario_cost_1950 = value_of_pension(fill(365, 40), reshaped_m2[61:100,:cooked_60_1950] ./ 100, 1.07, 40)
	scenario_cost_2019 = value_of_pension(fill(365, 40), reshaped_m2[61:100,:cooked_60_2019] ./ 100, 1.07, 40)
end

# ╔═╡ 2c6320da-fb5b-4bb4-ba8d-8f8f20a74e88
(scenario_cost_1950, scenario_var, scenario_cost_2019)

# ╔═╡ cdc7fa20-a8df-4d96-b201-c98a4ea7f6b4
md"""
The numbers you are seeing in these bar charts are the prices of **an asset that promises to pay 365 Rs per year to someone retiring at age 60**.

The asset is called an **annuity**.

So, if you went to an insurance company to buy an asset that paid you Rs365 a year, every year until you died, it would cost you $(round(scenario_cost_2019, digits=2)). In 1950, it would have cost you $(round(scenario_cost_1950, digits=2)). The asset you would be buying is called an annuity.

Insurance companies are in the business of pricing annuities. The price of the annuity depends on the insurance company's estimation that a typical person is going to live to an old-age.

As people live longer lives, the lifetime payout of the annuity has gone up. So it costs more to buy an annuity.

When the government promises a pension, it is essentially incurring a debt. It is taking upon itself a liability, whose cost *today* is just like the price of an annuity.

As people live longer, the cost of promising a pension keeps increasing.
"""

# ╔═╡ f97ce385-1e2d-4af0-9321-d88822de2fe4
md"""
But these are the costs per-person.

When the age structure of the population changes, and there are more old-age dependent people, the total pension burden increases.
"""

# ╔═╡ 7d5023cf-66f7-42d4-b754-536d21c37aae
# begin	
# 	# Create a slider with a default value and bind it to a variable
# 	@bind pension_amount_viz Slider(365:5000, default=365, show_value=true)
# end

# ╔═╡ 7e99e7dc-1294-4c2a-8d9a-e5f21492515c
md"""
## Age of retirement
"""

# ╔═╡ 911955e2-62a6-474f-bb5d-1c2fd4865b5b
md"""
Set the age of retirement
$(@bind retirement_age Slider(35:60, default=60, show_value=true))

Set the pension amount
$@bind pension_amount_viz Slider(365:5000, default=365, show_value=true)
"""

# ╔═╡ 8a7843a3-a719-4317-8999-7def2ebecc93
begin
	# minimum_pension = 365
	max_expected_lifespan = 100
	# retirement_age = 60
	years_after_retirement = max_expected_lifespan - retirement_age
	discount_rate = 1.07
end

# ╔═╡ 74e31918-8fc3-4ac5-bd8b-89965269aa70
p_mat=make_pension_mat(years_after_retirement, index, pension_amount);

# ╔═╡ c13efbff-b9c8-41a2-b31a-e0398dfd2fe7
retirement_age, years_after_retirement, index, pension_amount

# ╔═╡ 0fecad2b-9933-4cda-a3d2-2c64d234d1e0
# md"""
# # $(@bind c Button("🎖️"))
# # $(@bind indexation Select(["prices" => "🍞", "wages" => "🧑🏽‍🤝‍🧑🏽"]))
# """

# ╔═╡ cfde1e99-3b6b-4012-95d6-ce5ba79ebd85
md"""Lets calculate the present value of the pension for someone retiring at $retirement_age with a pension of Rs $pension_amount_viz a year.)"""
# The pension will be indexed to $indexation
# The maximum expected lifespan is $max_expected_lifespan""")

# ╔═╡ 18905212-226c-44c0-b880-6973fb7f217e
md"""
### Together
"""

# ╔═╡ 7edc861c-67a5-4b54-9011-e0bec6be78c1
plot(fig1, fig2, layout=(1, 2))

# ╔═╡ dfbffc17-4fc0-4e58-a9de-d4fd40159357
md"""
_insert the old-age dependency ration, with projections going up to 2100_
"""

# ╔═╡ 5989dbf6-35df-454f-99ce-e14b71590d48
md"""
There are a few different questions here:
- Retirement age
- Indexation: Price indexed versus wage-indexed annuities
- Debt burden: It is very difficult to get you head around the magnitude of the pension burden. So don't be tempted to increase the pension from 365 to 5000. It does not help in understanding how high a large number is. make it relative. What is the actual pension burden going to be on the state governments, as a fraction of their revenue income? 
- What are the projections regarding India's age structure?
"""

# ╔═╡ b5f353c5-772a-46df-b239-2e107b4f68b1


# ╔═╡ a1b91341-e477-4f78-a3d1-6938a12908e4


# ╔═╡ 46d98c8e-7898-466c-8549-ae8f6e65ebde


# ╔═╡ 5b9e1f99-2051-4927-a6c9-e369ff45d0ab


# ╔═╡ d2301db8-37a6-46f2-8f50-27b2d26a862f
begin
	col2 = Symbol("cooked_$retirement_age")
	value_of_pension2(p_mat, mortality_rate[retirement_age+1:100,col2] ./ 100, 1.07, years_after_retirement)
end

# ╔═╡ f37b1c40-0090-4efb-8966-730293a81322
# ╠═╡ disabled = true
#=╠═╡
begin
	col2 = Symbol("cooked_$retirement_age")
	value_of_pension2(p_mat, mortality_rate[retirement_age+1:100,col2] ./ 100, 1.07, years_after_retirement)
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CSV = "~0.10.10"
DataFrames = "~1.5.0"
PlotThemes = "~3.1.0"
Plots = "~1.38.13"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "01453e86dff6f8b6e9e8348958f1ea79b177f939"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "ed28c86cbde3dc3f53cf76643c2e9bc11d56acc7"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.10"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8b8a2fd4536ece6e554168c21860b6820a8a83db"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "19fad9cd9ae44847fe842558a744748084a722d1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e77dbf117412d4f164a464d610ee6050cc75272"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.6"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "099e356f267354f46ba65087981a77da23a279b7"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.0"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a5aef8d4a6e8d81f171b2bd4be5265b01384c74c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.10"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "3c5106dc6beba385fd1d37b9bf504271f8bfa916"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.13"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "213579618ec1f42dea7dd637a42785a608b1ea9c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "ba4aa36b2d5c98d6ed1f149da916b3ba46527b2b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.14.0"

    [deps.Unitful.extensions]
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═620e83fe-0e29-4217-85eb-7cb7284ef499
# ╠═64eb2201-2af5-4b06-bc76-cef60ba0cad4
# ╟─2fd5e8fe-c195-4a71-8f42-2260bad5ce0a
# ╟─9a00fe06-9828-4208-b655-58e933aebc08
# ╟─ab4bd33c-9e04-427d-8808-fdfc2ea56742
# ╟─3de95442-fca6-11ed-17dd-41fc8582cd5c
# ╟─b27f4c1d-8ee4-4d0b-97ab-380d89dacd28
# ╟─06384b92-5333-4b5f-849a-fb13ff6b571e
# ╠═2b76ad55-d09f-4fa8-9a98-ea3ff32b3ac6
# ╠═ebd0ce86-e18f-4d22-bcdf-d88c2b9d5b3f
# ╟─277797ab-724a-4848-afb3-1d2265ffde48
# ╠═feb4581a-c97e-4259-a460-531185d28326
# ╟─b331d2ba-67da-44cd-8101-f1ec6905d77f
# ╟─ae898e69-1ec1-4ca1-943c-99a8af3096cc
# ╠═f3579a67-2d9b-4e18-91b9-54aa806c4418
# ╟─ea710959-34dd-49f4-851b-a029c045c757
# ╟─06351924-c4d1-4ab8-bde3-3a3504677244
# ╠═0f5333ef-3039-4470-8886-190ee368dd92
# ╠═0d141036-2508-4fc5-a0fd-f23109638f3f
# ╠═f566dfb2-66da-40f3-bf03-0992d769873d
# ╟─8a7843a3-a719-4317-8999-7def2ebecc93
# ╟─2abd32aa-04b1-46e3-bf5b-740530773a2c
# ╟─60dbedcb-01f1-454b-9146-0a40785be747
# ╠═c13efbff-b9c8-41a2-b31a-e0398dfd2fe7
# ╠═74e31918-8fc3-4ac5-bd8b-89965269aa70
# ╠═6914def4-a379-4b3b-8ebe-d307b7d3ca1f
# ╠═f37b1c40-0090-4efb-8966-730293a81322
# ╠═5370a609-07d7-4ef5-8dbe-361e25fbe853
# ╠═cf4e0c5d-2bdb-4a1c-895d-e7d3570627c7
# ╠═4030a33b-5771-4fb5-b8f8-53496cfe13ad
# ╠═4815ae30-f791-4bb3-8543-8d96193c461e
# ╠═cb87055a-a240-4129-98f6-6d87a6bdd8fa
# ╠═e7d550b2-be57-41d9-89aa-ffa4980eae3f
# ╟─345ce785-fc32-4340-8d42-e557bfe5d3f1
# ╠═c9351453-1d94-4360-b393-4e74b0640ba9
# ╠═fb427276-0d24-4570-aced-f457c9477b02
# ╠═b4b4ca0d-8a1f-4ea2-8250-73015865c226
# ╠═82d54039-bca5-4c3b-8ef1-8ad88958c777
# ╠═4ae1dbcc-0938-4605-91e6-9c0da8a4c85b
# ╠═3704356c-44fd-4bd3-87f8-3652f83de810
# ╠═d5ec763f-6021-4566-a16f-4c219fe5116e
# ╠═f930becd-d683-434f-b97f-3fe81a1f8e8b
# ╠═185ae1c3-a68b-415d-b39c-be751c0779ff
# ╠═92d52b8e-3436-45b6-95ac-6c57602d3c30
# ╠═c40a18ce-5bfe-4ad0-91b0-b0454136abf3
# ╠═d6344906-be05-4f49-9f19-9b770cd36d35
# ╠═a077cc3b-d183-4528-b5ba-3c924bf929f1
# ╠═babc50fd-c355-4ea9-a749-eb206d0e3495
# ╠═072efd53-a9e5-4ad5-bef3-bb3bcb33b3eb
# ╠═4d22d4c3-ef7e-4e35-8f95-645830f8929b
# ╠═ba0b0292-8188-4797-8185-95ffc29b0d25
# ╠═222adafb-24bb-443c-be0c-bee6d69f641d
# ╠═a804f388-31de-43db-ba0e-600d02c26f06
# ╠═73a8b5ea-853a-41db-85f8-57ab8936996c
# ╠═c27ee5c1-997f-44a9-bc40-876bdc320e55
# ╠═0cb6f88d-7b8a-4d82-9795-310563aa31b4
# ╠═55da71d5-f553-4b05-9109-af83eadbaad6
# ╠═8a7bdd87-9a3a-404e-8a16-71ead8074072
# ╠═e28f31bd-d229-4e0e-a48f-6854493fedf5
# ╠═a7a17973-c45b-4826-b156-b20f4ea28b2f
# ╟─6f752523-f2e2-4321-9453-3de46d5b4169
# ╟─ae181d07-096c-4428-831a-5ee805810fef
# ╠═f41dcd9f-9611-4ed5-92ba-6c5009221fcb
# ╟─15e4be92-93a7-4273-9e32-856006ddcb61
# ╟─813f506d-1777-49a7-afd7-e962588b6290
# ╟─6be6c3c0-3974-42d6-a99e-213fd948bda5
# ╟─bc867c54-b9b4-4b72-97c8-85af8136afb0
# ╟─5213433c-9913-4d83-8b9e-83c09fbac19c
# ╠═51a13f1a-d305-4379-a23d-dc1affabeb92
# ╠═f040371c-f4a7-467c-90cb-d32ede5c6d2f
# ╟─2667661c-4806-4a3a-8e56-3691d2ac93ce
# ╠═2b046c22-5b07-4b31-8ec5-351817273a0c
# ╠═2c6320da-fb5b-4bb4-ba8d-8f8f20a74e88
# ╟─cdc7fa20-a8df-4d96-b201-c98a4ea7f6b4
# ╠═f97ce385-1e2d-4af0-9321-d88822de2fe4
# ╠═7d5023cf-66f7-42d4-b754-536d21c37aae
# ╟─7e99e7dc-1294-4c2a-8d9a-e5f21492515c
# ╠═d2301db8-37a6-46f2-8f50-27b2d26a862f
# ╟─911955e2-62a6-474f-bb5d-1c2fd4865b5b
# ╟─0fecad2b-9933-4cda-a3d2-2c64d234d1e0
# ╠═cfde1e99-3b6b-4012-95d6-ce5ba79ebd85
# ╟─18905212-226c-44c0-b880-6973fb7f217e
# ╠═7edc861c-67a5-4b54-9011-e0bec6be78c1
# ╟─dfbffc17-4fc0-4e58-a9de-d4fd40159357
# ╠═5989dbf6-35df-454f-99ce-e14b71590d48
# ╠═b5f353c5-772a-46df-b239-2e107b4f68b1
# ╠═a1b91341-e477-4f78-a3d1-6938a12908e4
# ╠═46d98c8e-7898-466c-8549-ae8f6e65ebde
# ╠═5b9e1f99-2051-4927-a6c9-e369ff45d0ab
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

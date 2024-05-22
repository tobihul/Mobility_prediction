using CSV, GLM, DataFrames, StatsPlots, TypedTables, Statistics,LinearAlgebra, Distributions

data = CSV.read("C:\\Users\\uqthulle\\Downloads\\housingdata.csv", DataFrame)

X = data[:,1]

Y = data[:,2]

data = Table(X = X, Y = Y)

scatter(X,Y, xlabel = "House price (\$)", ylabel = "House size (m²)", label = false)
# Perform linear regression
model = lm(@formula(Y ~ X), data)

y_hat = coef(model)[2].*X.+coef(model)[1]

x = [ones(length(X),1) X]
y = Y
b = pinv(x'*x)*x'*y # Regressed parameters
y_hat = x*b

# Now we plot the model
plot!(X,y_hat, label = false, c = :black)




############################# PART B ###########################

# We calculate the MS_res from our actual model.
df_m1 = size(x,1)-size(x,2) # This is the degrees of freedom. Note that it is always the number of rows of x minus the number of columns of x.
MS_res = sum((y.-y_hat).^2)./df_m1 # Mean squares due to regression

# Checking for the significance of b1.
# We fit a constant first, P=b0
x_m0 = ones(length(X),1)
b_m0 = pinv(x_m0'*x_m0)*x_m0'*y
y_m0_hat = x_m0*b_m0

df_m0 = size(x_m0,1)-size(x_m0,2) # This is the degrees of freedom. Note that it is always the number of rows of x minus the number of columns of x.
dv = df_m0 - df_m1
MS_addingB1 = sum((y_hat.-y_m0_hat).^2)/dv# Obviously, I could skip this "division by 1", but I put it ther for clarity (degrees of freedom = 1)


F = MS_addingB1/MS_res
pvalue = 1-cdf(FDist(dv, df_m1),F) # Fitting b1 is sensible at alpha 0.10.

############################# PART C ###########################

# Residuals' plot
plot(X,y-y_hat,seriestype = :scatter) # I don't see nothing strange there (no need to pass an official Grubbs' test)

# Error in b1
se2 = MS_res # Already claculated above
vc = se2*pinv(x'*x) # This is the variance-covariance matrix
b1 = b[2] # the second row of the b-vector is the b1 parameter (the first is the offset, b0).
s_b1 = sqrt(vc[2,2]) # The element 2,2 is the variance in b1. The square root transforms it into standard deviation

# 95 CI for b1
alpha_val = 0.05
t = quantile(TDist(df_m1),1-alpha_val/2)
CI = [b1-t*s_b1 b1+t*s_b1]

############################# PART D ###########################

x0 = [1 0.5]
y_hat0 = x0*b                               # This is the expected value at T=290
s_y_hat = sqrt(x0*vc*x0')
CI = [y_hat0-t*s_y_hat y_hat0+t*s_y_hat]    # this is the confidence interval for the y_hat

s_y = sqrt(se2.+x0*vc*x0')
CI = [y_hat0-t*s_y y_hat0+t*s_y]            # This is the confidence interval for y.

# BONUS: If you want, we can also create the full CI lines!
xi=collect(minimum(X):1:maximum(X))
y_hat0=zeros(length(xi),1)
s_y_hat=zeros(length(xi),1)
CI_left=zeros(length(xi),1)
CI_right=zeros(length(xi),1)
for i=1:length(xi)
    x0=[1 xi[i]]
    y_hat0[i] = (x0*b)[1]                   # This is the expected value at T=290
    s_y_hat[i] = sqrt(x0*vc*x0')[1]
    CI_left[i] = (y_hat0[i]-t*s_y_hat[i])[1]    # this is the confidence interval for the y_hat
    CI_right[i] = (y_hat0[i]+t*s_y_hat[i])[1]
end

plot!(xi, [CI_left CI_right], c = :black, label = false)

# BONUS: If you want, we can also create the full CI lines for the experimental response!
xi=collect(minimum(X):1:maximum(X))
s_y=zeros(length(xi),1)
CI_exp_left=zeros(length(xi),1)
CI_exp_right=zeros(length(xi),1)
for i=1:length(xi)
    x0=[1 xi[i]]
    s_y[i]=sqrt(se2.+x0*vc*x0')[1]
    CI_exp_left[i]=y_hat0[i]-t*s_y[i]  # this is the confidence interval for the y_hat
    CI_exp_right[i]=y_hat0[i]+t*s_y[i]
end
plot!(xi, [CI_exp_left, CI_exp_right], linestyle = :dash, c = 1, label = false)

function linear_regression(X,Y)
    x = [ones(length(X),1) X]
    y = Y
    b = pinv(x'*x)*x'*y # Regressed parameters
    y_hat = x*b 
    df_m1 = size(x,1)-size(x,2)
    MS_res = sum((y.-y_hat).^2)./df_m1
    se2 = MS_res # Already claculated above
    vc = se2*pinv(x'*x)
    alpha_val = 0.05
    t = quantile(TDist(df_m1),1-alpha_val/2)
    xi=collect(minimum(X):1:maximum(X))
    y_hat0=zeros(length(xi),1)
    s_y_hat=zeros(length(xi),1)
    CI_left=zeros(length(xi),1)
    CI_right=zeros(length(xi),1)
    for i=1:length(xi)
        x0=[1 xi[i]]
        y_hat0[i] = (x0*b)[1]                   # This is the expected value at T=290
        s_y_hat[i] = sqrt(x0*vc*x0')[1]
        CI_left[i] = (y_hat0[i]-t*s_y_hat[i])[1]    # this is the confidence interval for the y_hat
        CI_right[i] = (y_hat0[i]+t*s_y_hat[i])[1]
    end
    xi=collect(minimum(X):1:maximum(X))
    s_y=zeros(length(xi),1)
    CI_exp_left=zeros(length(xi),1)
    CI_exp_right=zeros(length(xi),1)
    for i=1:length(xi)
        x0=[1 xi[i]]
        s_y[i]=sqrt(se2.+x0*vc*x0')[1]
        CI_exp_left[i]=y_hat0[i]-t*s_y[i]  # this is the confidence interval for the y_hat
        CI_exp_right[i]=y_hat0[i]+t*s_y[i]
    end
    scatter(X,Y, xlabel = "X", ylabel = "Y", label = false)
    plot!(X, y_hat, label = false, c = :black)
    plot!(xi, [CI_left CI_right], c = :blue, label = false, fillrange=(CI_left), fillalpha = 0.1) # fill the area between the confidence lines
    p = plot!(xi, [CI_exp_left, CI_exp_right], linestyle = :dash, c = :black, label = false, fillrange=(CI_exp_left), fillalpha = 0.1)
        return p
end

scatter(X,Y, xlabel = "X", ylabel = "Y", label = false)
plot!(X, y_hat, label = false, c = :black)
plot!(xi, [CI_left CI_right], c = :blue, label = false, fillrange=(CI_left), fillalpha = 0.1) # fill the area between the confidence lines
p = plot!(xi, [CI_exp_left, CI_exp_right], linestyle = :dash, c = :black, label = false, fillrange=(CI_exp_left), fillalpha = 0.1)



cor(X,Y)


p = linear_regression(X,Y)



cor(X,Y)
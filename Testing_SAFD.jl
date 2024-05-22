using MS_Import, SAFD, LoopVectorization, StatsPlots

function mass_align(Rt::Vector{Float32}, Mz_values::Matrix{Float32}, Mz_intensity::Matrix{Float32})
    # Round the MS-values to the nearest integer
    Mz_values = round.(Mz_values .* 10) ./ 10
    # Get the unique masses (that were just rounded to the nearest integer)
    unique_mz_values::Vector{Float32} = sort(unique(Mz_values))
    # Define the intensity matrix
    plot_matrix::Matrix{Float32} = zeros(length(Rt), length(unique_mz_values))
    # Make variables to save memory allocations
    inner_mz_values = Vector{Float32}(undef, size(Mz_values, 2))
    inner_mz_intensity = Vector{Float32}(undef, size(Mz_values, 2))
    inner_plot_matrix = Vector{Float32}(undef, size(plot_matrix, 2))
    i = 1
    for i in 1:length(Rt)
        # Saving a view of the matrix (less memory allocation)
        inner_mz_values .= view(Mz_values, i, :)
        inner_mz_intensity .= view(Mz_intensity, i, :)
        inner_plot_matrix .= view(plot_matrix, i, :)
        # Set pos to 1, start of Julia counting in arrays
        pos = 1
        #last_pos = 1 # removed code
        # Make a turbo loop to increase speed
        @tturbo for k in 1:length(inner_mz_values)
            # Loop over the MZ values (k) for this retention moment (i)
            for j in 1:length(unique_mz_values)
                # Check for all the unique masses (j)
                #pos != last_pos && break
                # If the current mass is equal to the current unique mass, safe the intensity
                test_outcome = (unique_mz_values[j] == inner_mz_values[k])
                inner_plot_matrix[j] += test_outcome * inner_mz_intensity[k]
                pos = j * test_outcome + !test_outcome
            end
        end
        plot_matrix[i, :] .= inner_plot_matrix
    end
    return unique_mz_values, plot_matrix
end

pathin = "C:\\Users\\uqthulle\\Documents\\Test data" 
filenames = ["2710 4pm_01_inj01_32.mzXML"]
mz_thresh = [0, 0] #Sets threshold for mz values
int_thresh = 500 #Remove most of the noise

GC.gc()

mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
polarity, Rt = import_files_MS1(pathin, filenames, mz_thresh, int_thresh)
FileName = m[1]
#Adjust the settings for SAFD here
max_numb_iter = 1000 #Roughly the number of features that will be found, if there are less than n_iter it will stop automatically
max_t_peak_w=300 # or 20
res=20000
min_ms_w=0.02
r_thresh=0.9 #How well the fitted gaussian overlaps the original peak
min_int=1000
sig_inc_thresh=5
S2N=3 #minimum signal to noise ratio for SAFD to take data as a feature
min_peak_w_s=3

unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
GC.gc()
#Run SAFD (this may take some time depending on n iterations and min_int)
rep_table, SAFD_output = safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)

heatmap(Rt, unique_mz_values, plot_matrix',
        #c = cgrad([:white,:navy,:indigo,:teal,:green,:yellow,:red],[0,0.04,1]),
        color=:plasma,
        clims=(20000, 225000),
        size=(1280, 720),
        xlabel="Rt (min)",
        ylabel="m/z",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm,
        colorbar = false,
        dpi = 300)

        p = scatter!(SAFD_output[:, 4], SAFD_output[:, 8],
        #series_annotations = text.(1:length(SAFD_output[:,4]),size = 1),
        markershape=:xcross,
        c = :green,
        legend=:topleft,
        markersize=2.5,
        #title="$(filenames[1]), $max_numb_iter iterations, S/N = $S2N, r = $r_thresh, accepted_res = 1.5, Componetization -> ($(length(SAFD_output[:,1])) features)",
        left_margin=5Plots.mm, right_margin=17.5Plots.mm,
        bottom_margin=8.5Plots.mm,dpi = 300
    )

    surface(Rt, unique_mz_values, plot_matrix',
    c =:blues,
    clims=(0, 25000000))


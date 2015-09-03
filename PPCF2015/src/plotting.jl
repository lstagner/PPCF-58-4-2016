function make_stacked_plot{T<:Real}(energy::Array{T,1},pitch::Array{T,1},x::Array{T,2},err::Array{T,2},dist::Array{T,2}; cb_range = (), cmap = "gist_earth", filename = Void)

  # nr -> Pitch, nc -> Energy
  nr,nc = size(x)
  inds = 2:2:nc
  nrow= length(inds)
  cm = plt.cm[:get_cmap](cmap)
  cm[:set_under]("w")

  fig, ax = plt.subplots(nrows=nrow,sharex=true)
  fig[:set_size_inches]((5,5))
  if length(cb_range) != 2
    v_min,v_max = extrema(x)
  else
    if any(collect(cb_range) .== Void)
      if cb_range[1] == Void
        v_min,v_max = minimum(x), cb_range[2]
      else
        v_min,v_max = cb_range[1], maximum(x)
      end
    else
      v_min,v_max = cb_range
    end
  end
  p=0
  for (i,ind) in enumerate(inds)
    ax[nrow+1 - i][:set_frame_on](false)
    ax[nrow+1 - i][:get_xaxis]()[:set_ticks_position]("none")
    ax[nrow+1 - i][:fill_between](pitch,x[:,ind] .- err[:,ind] ,x[:,ind] .+ err[:,ind],facecolor="lightgray",color="lightgray")
    p = ax[nrow+1 - i][:scatter](pitch,dist[:,ind],c = dist[:,ind],vmin=v_min,vmax=v_max,cmap=cm,marker="o",linewidth=0.5)
    ax[nrow+1 - i][:set_xticks]([-1.0,-0.5,0.0,0.5,1.0])
    ax[nrow+1 - i][:set_yticks]([])
    #ax[nrow+1 - i][:set_ylabel](@sprintf("%d keV",energy[ind]),rotation=0.0,fontsize="small")
    ax[nrow+1 - i][:text](-1.0,dist[1,ind],@sprintf("%d keV -",energy[ind]),fontsize="small", va="center",ha="right")
  end

  ax[end][:set_xlabel]("Pitch",fontsize="small")
  ax[end][:tick_params](labelsize="small")

  fig[:tight_layout](h_pad=0.0)
  cax,kw = matplotlib[:colorbar][:make_axes](ax,aspect=40,shrink=0.95,pad=0.0)
  cb = fig[:colorbar](p, cax=cax,ticks=[])
  cb[:solids][:set_edgecolor]("face")

  if filename != Void
    fig[:savefig](filename,bbox_inches="tight",dpi=1200)
  end

  gc()
end

function make_ep_plot{T<:Real}(energy::Union(Array{T,1},Array{T,2}),pitch::Union(Array{T,1},Array{T,2}),x::Array{T,2}; cmap = "gist_earth",cb_range = (), filename = Void)

  # nr -> Pitch, nc -> Energy
  nr,nc = size(x)

  cm = plt.cm[:get_cmap](cmap)
  cm[:set_under]("w")

  if length(cb_range) != 2
    v_min,v_max = extrema(x)
  else
    if any(collect(cb_range) .== Void)
      if cb_range[1] == Void
        v_min,v_max = minimum(x), cb_range[2]
      else
        v_min,v_max = cb_range[1], maximum(x)
      end
    else
      v_min,v_max = cb_range
    end
  end

  fig,ax = plt.subplots()
  fig[:set_size_inches]((5,3.5))
  p=ax[:pcolor](energy,pitch,x,cmap = cm)
  ax[:set_xlim]((0.0,100.0))
  ax[:set_yticks]([-1.0,-0.5,0.0,0.5,1.0])
  ax[:set_xticks]([0.0,20.0,40.0,60.0,80.0,100.0])
  ax[:tick_params](labelsize="small")
  ax[:set_xlabel]("Energy [keV]",fontsize="small")
  ax[:set_ylabel]("Pitch",fontsize="small",labelpad=0)
  cb = fig[:colorbar](p,pad=0)
  cb[:ax][:tick_params](labelsize="small")
  cb[:ax][:yaxis][:get_offset_text]()[:set_size]("small")
  cb[:ax][:yaxis][:get_offset_text]()[:set_position]((-1.2,1))
  cb[:solids][:set_edgecolor]("face")
  fig[:tight_layout]()

  if filename != Void
    fig[:savefig](filename,bbox_inches="tight",dpi=1200)
  end

  gc()
end

function t_colormask(clip c, int color, int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	c.t_24colormask(color, tolerance=tolerance, expandblur=expandblur, BT601=BT601, blur=blur)
}

function t_2colormask(clip c, int color1, int "color2", int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	c.t_24colormask(color1, color2, tolerance=tolerance, expandblur=expandblur, BT601=BT601, blur=blur)
}

function t_3colormask(clip c, int color1, int "color2", int "color3", int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	c.t_24colormask(color1, color2, color3, tolerance=tolerance, expandblur=expandblur, BT601=BT601, blur=blur)
}

function t_6colormask(clip c, int color1, int "color2", int "color3", int "color4", int "color5", int "color6", int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	c.t_24colormask(color1, color2, color3, color4, color5, color6, tolerance=tolerance, expandblur=expandblur, BT601=BT601, blur=blur)
}

function t_9colormask(clip c, int color1, int "color2", int "color3", int "color4", int "color5", int "color6", int "color7", int "color8", int "color9", int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	c.t_24colormask(color1, color2, color3, color4, color5, color6, color7, color8, color9, tolerance=tolerance, expandblur=expandblur, BT601=BT601, blur=blur)
}

function t_12colormask(clip c, int color1, int "color2", int "color3", int "color4", int "color5", int "color6", int "color7", int "color8", int "color9", int "color10", int "color11", int "color12", int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	c.t_24colormask(color1, color2, color3, color4, color5, color6, color7, color8, color9, color10, color11, color12, tolerance=tolerance, expandblur=expandblur, BT601=BT601, blur=blur)
}

function t_24colormask(clip c, int color1, int "color2", int "color3", int "color4", int "color5", int "color6", int "color7", int "color8", int "color9", int "color10", int "color11", int "color12", int "color12", int "color13", int "color14", int "color15", int "color16", int "color17", int "color18", int "color19", int "color20", int "color21", int "color22", int "color23", int "color24", int "tolerance", int "expandblur", bool "BT601", int "blur")
{
	color2=default(color2,-1)
	color3=default(color3,-1)
	color4=default(color4,-1)
	color5=default(color5,-1)
	color6=default(color6,-1)
	color7=default(color7,-1)
	color8=default(color8,-1)
	color9=default(color9,-1)
	color10=default(color10,-1)
	color11=default(color11,-1)
	color12=default(color12,-1)
	color13=default(color13,-1)
	color14=default(color14,-1)
	color15=default(color15,-1)
	color16=default(color16,-1)
	color17=default(color17,-1)
	color18=default(color18,-1)
	color19=default(color19,-1)
	color20=default(color20,-1)
	color21=default(color21,-1)
	color22=default(color22,-1)
	color23=default(color23,-1)
	color24=default(color24,-1)
	
    tolerance=default(tolerance,10)
	expandblur=default(expandblur,0)
	BT601=default(BT601,false)
	blur=default(blur,1)

    colors_string = string(color1) +
            \ (color2 == -1 ? "" : " " + string(color2)) +
            \ (color3 == -1 ? "" : " " + string(color3)) +
            \ (color4 == -1 ? "" : " " + string(color4)) +
            \ (color5 == -1 ? "" : " " + string(color5)) +
            \ (color6 == -1 ? "" : " " + string(color6)) +
            \ (color7 == -1 ? "" : " " + string(color7)) +
            \ (color8 == -1 ? "" : " " + string(color8)) +
            \ (color9 == -1 ? "" : " " + string(color9)) +
            \ (color10 == -1 ? "" : " " + string(color10)) +
            \ (color11 == -1 ? "" : " " + string(color11)) +
            \ (color12 == -1 ? "" : " " + string(color12)) +
            \ (color13 == -1 ? "" : " " + string(color13)) +
            \ (color14 == -1 ? "" : " " + string(color14)) +
            \ (color15 == -1 ? "" : " " + string(color15)) +
            \ (color16 == -1 ? "" : " " + string(color16)) +
            \ (color17 == -1 ? "" : " " + string(color17)) +
            \ (color18 == -1 ? "" : " " + string(color18)) +
            \ (color19 == -1 ? "" : " " + string(color19)) +
            \ (color20 == -1 ? "" : " " + string(color20)) +
            \ (color21 == -1 ? "" : " " + string(color21)) +
            \ (color22 == -1 ? "" : " " + string(color22)) +
            \ (color23 == -1 ? "" : " " + string(color23)) +
            \ (color24 == -1 ? "" : " " + string(color24))
            
    c.tcolormask(colors_string, tolerance, BT601)
    expandblur > 0 ? binomialblur(expandblur,u=1,v=1).mt_lut("x 5 *") : last
    blur == 1 ? removegrain(20, -1) :  \
        blur>0 ? binomialblur(blur,u=1,v=1) : last
}

function t_colormask2(clip c, string colors, int "tolerance", bool "BT601", int "blur", int "expandblur", bool "gray", int "lutthr", bool "mt")
{
    expandblur = default(expandblur, 0)
    blur = default(blur, 1)
    gray = default(gray, false)
    
    c.tcolormask(colors, tolerance=tolerance, bt601=bt601, lutthr=lutthr, mt=mt )
    expandblur > 0 ? binomialblur(expandblur,u=1,v=1).mt_lut("x 5 *") : last
    blur == 1 ? removegrain(20, -1) :  \
        blur > 0 ? binomialblur(blur,u=1,v=1) : last
    gray ? mt_lut(chroma="128") : last
}
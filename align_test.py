from SimpleCV import *
import cv2
import cv
import numpy

def stitch(img1, img2):
    #create an empty destination image with combined size
    h1,w1 = img1.size()
    h2,w2 = img2.size()
    dst_img = Image((h1+h2,w1+w2))
    dst_array = numpy.array(dst_img.getMatrix())

    ster = StereoImage(img1,img2)
    H = ster.findHomography()[0]
    #homo = numpy.rot90(numpy.matrix(H))
    #transform to one image
    img2_array = numpy.array(img2.getMatrix())
    res_array = cv2.warpPerspective(src   = img2_array,
                                    M     = H,
                                    dsize = dst_img.size(),
                                    dst   = dst_array,
                                    flags = cv2.INTER_CUBIC,
                                   )
                                   
    #res_img = Image(res_array,colorSpace = ColorSpace.RGB).toBGR()
    res_img = Image(res_array)
    res_img = res_img.rotate90()
    # blit the img1 now on coordinate (0, 0).
    res_img = res_img.blit(img1, alpha=0.4)
    return res_img
    
################################################################################
# TEST CODE
################################################################################    
if __name__ == "__main__":

    IM1 = "test_images/a0.png"
    IM2 = "test_images/a1.png"
    #IM1 = "test_images/wave1.jpg"
    #IM2 = "test_images/wave2.jpg"

    img1 = Image(IM1)
    img2 = Image(IM2)
    res  = stitch(img1,img2)
    res.save("stitched.png")
